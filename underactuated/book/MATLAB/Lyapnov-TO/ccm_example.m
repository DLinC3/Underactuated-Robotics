function ccm_example
    clear all;
    close all;
    
    use_lqr = false;
    
    xd=[0;0;0];%desired final state
    [A0, B0] = gradients(xd);%linearized dynamics about xd
    Q= eye(3);% state cost matrix
    R = 1;% input cost matrix
    [K, S] = lqr(A0,B0,Q,R);%generate lqr controller for linearized dynamics
    
    % x0=[0.5;0.5;0.5];%initial conditions
    x0=[9;9;9];%initial conditions
    % x0=[6;0.5;6];%initial conditions
    
    x = sdpvar(1,3); %state variables
    y = sdpvar(3,1); %variables for ccm 
    c = sdpvar(1,9*6+3);% variables for dual contraction metric
    %matrices for ccm 
    W0 = reshape(c(1:9),3,3);
    W1 = reshape(c(10:18),3,3);
    W2 = reshape(c(19:27),3,3);
    W3 = reshape(c(28:36),3,3);
    W4 = reshape(c(37:45),3,3);
    W5 = reshape(c(46:54),3,3);
    rho = c(55:57);
    
    %nonlinear polynomial dynamics
    f = dynamics(x);
    
    % polynomail gradients 
    [A, B] = gradients(x);
    
    %construct dual metric (fill in)
    W = W0 + W1*x(1) + W2*x(1)^2 + W3*x(2) + W4*x(2)^2 + W5*x(1)*x(2);


    %construct multiplier (fill in)
    rho0 = evalRho(rho, x);
    
    %compute Wdot (fill in)
    Wdot = (W1+2*W2*x(1)+W5*x(2))*f(1,1)+(W3+2*W4*x(2)+W5*x(1))*f(2,1); 
    
    %contraction rate
    lambda = 0.5;
    
    % contraction condition (fill in)
    P = -Wdot + A*W + W*A' - rho0 * (B*B') + 2*lambda*W;
    
  
    % setup sos constraints and solve 
    %%%%%%%%%%%%%%% fill in here %%%%%%%%%%%%%%%%

    p_constraint = sos(y' * (-P) * y);
    

    eps_val = 1e-4;
    w_constraint = sos(y' * (W - eps_val*eye(3)) * y);
    

    rho_constraint = sos(rho0 - eps_val);
    

    constraints = [p_constraint, w_constraint, rho_constraint];
    

    solvesos(constraints, [], [], c);
    
    c_star = value(c);
        
    %extract values for rho (fill in)
    rhohat = c_star(55:57);
    
    %extract values for ccm matrices (fill in)
    What{1} = reshape(c_star(1:9),3,3);
    What{2} = reshape(c_star(10:18),3,3);
    What{3} = reshape(c_star(19:27),3,3);
    What{4} = reshape(c_star(28:36),3,3);
    What{5} = reshape(c_star(37:45),3,3);
    What{6} = reshape(c_star(46:54),3,3);
    
    T=10;% time horizon for the simulation
    dt =0.01; % time step
    N=floor(T/dt);% number of time steps
    
    x=zeros(3,N);% state variable
    x(:,1) = x0;% set initial condition
    xd = [0;0;0];
    for k=1:N
        k
        [f, g] =dynamics(x(:,k));% compute the nonlinear dynamics
        if(use_lqr)
            u = -K*x(:,k);%use lqr to compute the control input
        else
            u = CCMController(x(:,k),xd,What,rhohat,B); %use the CCM controller
        end
        xdot = f + g*u;% nonlinear dynamics
        x(:,k+1) = x(:,k) + xdot*dt;%integrate nonlinear dynamics
        xdis = x(:,k+1);
    end
    
    figure(1)
    hold on
    plot(x(1,:),'r')
    plot(x(2,:),'g')
    plot(x(3,:),'b')
    hold off

end

function u = CCMController(xi,xd,W0,rho0,B)
    [xs, dxds, ds]= computeGeodesic(xd,xi,W0);
    u = integrateDeltaK(xs,W0,rho0,B,dxds,ds);
end

%==========================================================================
% gradients for the polynomial dynamics (fill in)
%==========================================================================
function [dfdx, dfdu] = gradients(x)
dfdx = [-1, 0, 1;
        2*x(1) - 2*x(3), -1, -2*x(1) + 1;
        0, -1, 0];
dfdu = [0; 0; 1];
end

%==========================================================================
% nonlinear polynomial dynamics
%==========================================================================
function [f, g] = dynamics(x)
    f(1,1) = -x(1) + x(3);
    f(2,1)  = x(1)^2-x(2)-2*x(1)*x(3) + x(3);
    f(3,1) = -x(2);
    
    g = [0;0;1];
end

%==========================================================================
% integrate deltaK (fill in)
%==========================================================================
function K = integrateDeltaK(xs, W0, rho0, B, dxds, ds)  
    N = size(xs,2);
    K = 0;    
    for i = 1:N
        x_i = xs(:,i);
        W = evalDualMetric(W0, x_i);
        M = inv(W); 
        rho = evalRho(rho0, x_i);
        dK = -0.5 * rho * B' * M; 
        delta_K = dK * dxds(:,i) * ds;
        K = K + delta_K;
    end
end



%==========================================================================
% evaluate delta K (fill in)
%==========================================================================
function dK = evalDeltaK(rho,B,W)
dK = -0.5 * rho * B' * inv(W);
end


%==========================================================================
% eval rho (fill in)
%==========================================================================
function rho=evalRho(rho0,x)
    rho=rho0(1)+ rho0(2)*x(1)+rho0(3)*x(1)^2;
end

%==========================================================================
% evaluate the dual metric (fill in)
%==========================================================================
function Wx = evalDualMetric(W,x)
    Wx=W{1} + W{2}*x(1) + W{3}*x(1)^2 + W{4}*x(2) + W{5}*x(2)^2 + W{6}*(x(1)*x(2));
end

%==========================================================================
% evaluate the Metric from the dual Metric
%==========================================================================
function Mx = evalMetric(W0,x)
    Wx = evalDualMetric(W0,x);
    Mx = inv(Wx);
end

%==========================================================================
% cost function for computing the geodesic
%==========================================================================
function e = cost_function(X,W0)

N = (length(X)-1)/6;
x = reshape(X(1:N*3),3,N);
dxds = reshape(X(N*3+1:end-1),3,N);

%implement cost for computing geodesic

    
e = 0;  
for i = 1:N
    xi = x(:,i);
    Mx = evalMetric(W0, xi); 
    energy_term = dxds(:,i)' * Mx * dxds(:,i);
    e = e + energy_term;
end
e = e * 1 / N; % s at [0,1]


end

%==========================================================================
% function for computing the geodesic
%==========================================================================
function [xs, dxds, ds]= computeGeodesic(xi,xd,W0)
    N = 10;
    x0=zeros(3,N);
    m = (xd-xi)/N;
    x0(:,1) = xi;
    for i=2:N
        x0(:,i) = xi+i*m;%initialize geodesic
    end
    dxds=zeros(3,N);
    for i=1:N
        dxds(:,i) = [0.1;0.1;0.1];%initial derivative of curve
    end
    
    %bounds on the geodesic curve
    lb=-ones(3,N)*100;
    ub=ones(3,N)*100;
    
    % fill in the boundary conditions on the path x
    e = 0.25;
    lb(:,1) = xi - abs(xi)*e;
    ub(:,1) = xi + abs(xi)*e;
    lb(:,N) = xd - abs(xd)*e;
    ub(:,N) = xd + abs(xd)*e;
    
    %bounds of the gradient of the geodesic
    dxdslb=-ones(3,N)*100;
    dxdsub=ones(3,N)*100;
    
    x0 = [x0(:);dxds(:);1];
    lb = [lb(:);dxdslb(:);0];
    ub = [ub(:);dxdsub(:);100];
    options = optimset('Largescale','off','MaxFunEvals',30000);
    
    [x, fval] = fmincon(@(x)cost_function(x,W0), x0, [], [], [], [], lb, ub,@(x)constraint(x),options);
    
    xs = reshape(x(1:3*N),3,N);
    
    dxds = reshape(x(3*N+1:end-1),3,N);
    
    ds = x(end);
end

%==========================================================================
% direct transcription integration constraint
%==========================================================================
function [c,ceq] = constraint(X)
    N = (length(X)-1)/6;
    x = reshape(X(1:N*3),3,N);
    dxds = reshape(X(N*3+1:end-1),3,N);
    ds = X(end);
    ceq =[];
    c = [];
        for i=1:N-1
            x0 = x(:,i);
            x1 = x(:,i+1);
            dxds0 = dxds(:,i);
            dxds1 = dxds(:,i+1);
            dxdsc = (dxds0 + dxds1) / 2.0;
            ceq = [ceq;x0 - x1 + (ds / 6.0) * (dxds0 + 4 * dxdsc + dxds1)];
        end
end

