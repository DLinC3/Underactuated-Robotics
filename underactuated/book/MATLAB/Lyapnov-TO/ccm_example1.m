function ccm_example
    clear all;
    close all;
    
    use_lqr = false;
    
    % --- Linearization about the origin for LQR ---
    xd = [0;0;0]; % desired final state
    [A0, B0] = gradients(xd); % linearized dynamics about xd
    Q = eye(3);    % state cost matrix
    R = 1;         % input cost matrix
    [K, S] = lqr(A0,B0,Q,R); % generate LQR controller for linearized dynamics
    
    x0 = [9;9;9];  % initial conditions for simulation
    
    % --- Declare YALMIP variables for the dual metric parameters ---
    % There are 6 symmetric 3x3 matrices and a polynomial (degree 2 in x1)
    c = sdpvar(1, 9*6+3); % 6*9 = 54 + 3 = 57 decision variables
    W0 = reshape(c(1:9),3,3);
    W1 = reshape(c(10:18),3,3);
    W2 = reshape(c(19:27),3,3);
    W3 = reshape(c(28:36),3,3);
    W4 = reshape(c(37:45),3,3);
    W5 = reshape(c(46:54),3,3);
    % rho is a 3-element vector so that rho(x) = rho(1) + rho(2)*x1 + rho(3)*x1^2.
    rho = c(55:57);
    
    % --- Define state as a symbolic (column) variable ---
    x = sdpvar(3,1); % state variable (column vector)
    
    % --- Define the nonlinear dynamics and its gradients ---
    f = dynamics(x);
    [A, B] = gradients(x);
    
    % --- Construct the dual metric W(x) as a polynomial in x1 and x2 ---
    % (Note: x3 does not appear by design.)
    W = W0 + W1*x(1) + W2*(x(1)^2) + W3*x(2) + W4*(x(2)^2) + W5*x(1)*x(2);
    
    % --- Construct the multiplier: here we simply use the polynomial rho ---
    rho0 = rho(1) + rho(2)*x(1) + rho(3)*x(1)^2;
    
    % --- Compute Wdot = (dW/dx1)*f1 + (dW/dx2)*f2 ---
    % f1 = -x(1) + x(3) and f2 = x(1)^2 - x(2) - 2*x(1)*x(3) + x(3)
    Wdot = (W1 + 2*W2*x(1) + W5*x(2)) * (-x(1) + x(3)) + ...
           (W3 + 2*W4*x(2) + W5*x(1)) * (x(1)^2 - x(2) - 2*x(1)*x(3) + x(3));
       
    % --- Define the contraction rate ---
    lambda = 0.5;
    
    % --- Formulate the contraction LMI ---
    % The contraction condition is:
    % -Wdot + A*W + W*A' - rho(x)*B*B' + 2*lambda*W <= 0.
    % Note that we use the evaluated polynomial rho(x) via evalRho.
    P = -Wdot + A*W + W*A' - rho0* (B*B') + 2*lambda*W;
    
    % --- Set up SOS constraints ---
    sos_constr = [];
    v = sdpvar(3,1); % auxiliary vector for scalarization
    % Require that v'*P*v is SOS (i.e. P is negative semidefinite in a polynomial sense)
    sos_constr = [sos_constr, sos(v'*P*v)];
    % Also require that the metric is uniformly positive definite.
    eps_val = 1e-3;
    sos_constr = [sos_constr, sos(v'*(W - eps_val*eye(3))*v)];
    
    % --- Solve the SOS program using MOSEK ---
    options = sdpsettings('solver','mosek','verbose',1);
    sol = solvesos(sos_constr,[],options,c);
    
    % --- Extract the optimized dual metric coefficients and multiplier ---
    rhohat = value(rho);
    What{1} = value(W0);
    What{2} = value(W1);
    What{3} = value(W2);
    What{4} = value(W3);
    What{5} = value(W4);
    What{6} = value(W5);
    
    % --- Simulation parameters ---
    T = 10;       % time horizon for the simulation
    dt = 0.01;    % time step
    N_sim = floor(T/dt); % number of time steps
    x_sim = zeros(3, N_sim+1); % state trajectory for simulation
    x_sim(:,1) = x0; % initial condition
    xd = [0;0;0];    % target state
    
    % --- Simulate the nonlinear system with the CCM controller ---
    for k = 1:N_sim
        [f_val, g_val] = dynamics(x_sim(:,k));
        if use_lqr
            u = -K*x_sim(:,k);
        else
            u = CCMController(x_sim(:,k), xd, What, rhohat, B);
        end
        xdot = f_val + g_val*u;
        x_sim(:,k+1) = x_sim(:,k) + xdot*dt;
    end
    
    % --- Plot the state trajectories ---
    figure(1)
    hold on
    plot(x_sim(1,:),'r')
    plot(x_sim(2,:),'g')
    plot(x_sim(3,:),'b')
    hold off
    title('State trajectories')
    xlabel('Time step')
    ylabel('State values')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CCM Controller: computes control by integrating the differential feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = CCMController(xi, xd, Wcoeff, rho0, B)
    [xs, dxds, ds] = computeGeodesic(xd, xi, Wcoeff);
    u = integrateDeltaK(xs, Wcoeff, rho0, B, dxds, ds);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradients for the polynomial dynamics: compute Jacobian A(x) and B(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, B] = gradients(x)
    % x is assumed to be a column vector [x1;x2;x3].
    % f(x) = [ -x1 + x3;
    %          x1^2 - x2 - 2*x1*x3 + x3;
    %          -x2 ]
    % Therefore, we compute:
    A = [ -1,            0,         1;
          2*x(1) - 2*x(3), -1,   -2*x(1) + 1;
          0,            -1,         0 ];
    % From the dynamics we have g = [0; 0; 1]
    B = [0;0;1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear polynomial dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, g] = dynamics(x)
    % x can be symbolic (sdpvar) or numerical.
    f = zeros(3,1);
    f(1) = -x(1) + x(3);
    f(2) = x(1)^2 - x(2) - 2*x(1)*x(3) + x(3);
    f(3) = -x(2);
    
    g = [0;0;1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integrateDeltaK: integrates the differential control along the geodesic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = integrateDeltaK(xs, Wcoeff, rho0, B, dxds, ds)
    % xs: 3xN points along the geodesic
    % dxds: 3xN derivative of the geodesic w.r.t. the path parameter s
    % ds: total geodesic length (scalar)
    N = size(xs,2);
    delta_u = zeros(1,N);
    for i = 1:N
        x_i = xs(:,i);
        % Evaluate the dual metric at x_i using the cell array of coefficients.
        W_x = evalDualMetric(Wcoeff, x_i);
        W33 = W_x(3,3);
        % Evaluate the multiplier rho(x) at x_i.
        rho_val = evalRho(rho0, x_i);
        % Differential control: delta u = - (rho)/(2*W33)* (change in x3)
        delta_u(i) = - (rho_val)/(2*W33) * dxds(3,i);
    end
    % Integrate using the trapezoidal rule along s from 0 to ds.
    s_grid = linspace(0, ds, N);
    u = trapz(s_grid, delta_u);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evalDeltaK: evaluates the differential control law at a given state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dK = evalDeltaK(rho_val, B, W_x)
    % For a given state, the differential control law is
    % dK = - (rho_val)/(2*W33) * B', with B = [0;0;1]
    W33 = W_x(3,3);
    dK = - (rho_val)/(2*W33) * (B'); % results in a 1x3 row vector
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evalRho: evaluates the polynomial multiplier rho(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = evalRho(rho0, x)
    % rho(x) = rho0(1) + rho0(2)*x1 + rho0(3)*x1^2.
    r = rho0(1) + rho0(2)*x(1) + rho0(3)*(x(1)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evalDualMetric: evaluates the dual metric W(x) at a given state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Wx = evalDualMetric(Wcoeff, x)
    % Wcoeff is assumed to be a cell array {W0, W1, ..., W5}.
    % Our polynomial is: W(x) = W0 + W1*x1 + W2*x1^2 + W3*x2 + W4*x2^2 + W5*x1*x2.
    Wx = Wcoeff{1} + Wcoeff{2}*x(1) + Wcoeff{3}*(x(1)^2) + ...
         Wcoeff{4}*x(2) + Wcoeff{5}*(x(2)^2) + Wcoeff{6}*(x(1)*x(2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evalMetric: computes the Riemannian metric M(x) as the inverse of the dual metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mx = evalMetric(Wcoeff, x)
    Wx = evalDualMetric(Wcoeff, x);
    Mx = inv(Wx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cost_function: cost for computing the geodesic via direct transcription
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = cost_function(X, W0)
    % X contains both the state along the path and its derivative plus ds:
    % Let N be the number of discretization points.
    % First 3*N entries are states, next 3*N are derivatives, last entry is ds.
    N = (length(X)-1)/6;
    x = reshape(X(1:3*N),3,N);
    dxds = reshape(X(3*N+1:end-1),3,N);
    ds = X(end);
    e = 0;
    for i = 1:N
        Mx = evalMetric(W0, x(:,i));
        e = e + dxds(:,i)' * Mx * dxds(:,i);
    end
    e = e / N; %
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeGeodesic: computes the geodesic between the current state xi and target xd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xs, dxds, ds] = computeGeodesic(xi, xd, Wcoeff)
    N = 10;  % number of discretization points
    % Initialize the state path guess: linearly interpolate from xi to xd.
    x0 = zeros(3,N);
    m = (xd - xi) / (N-1);
    for i = 1:N
        x0(:,i) = xi + (i-1)*m;
    end
    % Initialize derivative guess (can be a small constant vector).
    dxds0 = repmat([0.1;0.1;0.1], 1, N);
    
    % --- Set bounds for the state trajectory ---
    lb_state = -ones(3,N)*100;
    ub_state = ones(3,N)*100;
    % Impose boundary conditions: fix first and last states.
    lb_state(:,1) = xi-0.3;
    ub_state(:,1) = xi+0.3;
    lb_state(:,N) = xd-0.3;
    ub_state(:,N) = xd+0.3;
    
    % Bounds for the derivative (loose bounds)
    lb_dxds = -ones(3,N)*100;
    ub_dxds = ones(3,N)*100;
    
    % Form the decision vector: states, derivatives, and ds (scalar).
    X0 = [x0(:); dxds0(:); 1];
    lb = [lb_state(:); lb_dxds(:); 0];
    ub = [ub_state(:); ub_dxds(:); 100];
    
    options = optimset('Largescale','off','MaxFunEvals',30000,'Display','iter');
    [X_opt, fval] = fmincon(@(X) cost_function(X, Wcoeff), X0, [], [], [], [], lb, ub, @(X) constraint(X), options);
    
    xs = reshape(X_opt(1:3*N), 3, N);
    dxds = reshape(X_opt(3*N+1:end-1), 3, N);
    ds = X_opt(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constraint: dynamic consistency constraint for the geodesic transcription
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, ceq] = constraint(X)
    N = (length(X)-1)/6;
    x_path = reshape(X(1:3*N), 3, N);
    dxds = reshape(X(3*N+1:end-1), 3, N);
    ds = X(end);
    ceq = [];
    c = [];
    % Use Simpson's rule to enforce the integration constraint along the path.
    for i = 1:N-1
        x0 = x_path(:,i);
        x1 = x_path(:,i+1);
        dxds0 = dxds(:,i);
        dxds1 = dxds(:,i+1);
        dxdsc = (dxds0 + dxds1) / 2.0;
        ceq = [ceq; x0 - x1 + (ds/6.0)*(dxds0 + 4*dxdsc + dxds1)];
    end
end
