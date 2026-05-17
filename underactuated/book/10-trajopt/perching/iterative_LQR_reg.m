%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% notes:
% iLQR with a trust-region-style quadratic shaping term:
%   - Maintains a nominal trajectory (xtraj0, utraj0).
%   - Backward pass computes local gains in deviations:
%         δu_k = k_k + K_k δx_k,    δx_k := x_k - x_k^{nom}
%   - Forward pass rolls out an updated nominal under the improved policy:
%         u_k^{new} = u_k^{nom} + α k_k + K_k (x_k^{new} - x_k^{nom})
%         x_{k+1}^{new} = x_k^{new} + f(x_k^{new}, u_k^{new}) dt
%   - The quadratic shaping term in cost() is centered at (x0,u0) = (xtraj00,utraj00),
%     i.e., it penalizes deviation from the previous nominal (a “trust region” proxy):
%         ℓ_k = 0.5 (x_k - x̄_k)^T Q (x_k - x̄_k) + 0.5 (u_k - ū_k)^T R (u_k - ū_k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
% iterativeLQR function
%==========================================================================
function [xtraj, utraj, ktraj, Ktraj] = iterative_LQR_reg( ...
    x0, xtraj, utraj, ktraj, Ktraj, N, dt, param, Q, R, Qf, xd, DYNAMICS, GRADIENTS)

    % J0 is treated as the previous cost for line-search acceptance
    J    = 1e6;
    Jlast = J;

    % Save the initial nominal as the trust-region anchor (x̄_k, ū_k)
    xtraj0 = xtraj;
    utraj0 = utraj;

    for i=1:1000
        % Forward pass: rollout an improved nominal using (k,K) and step size α
        [xtraj, utraj, J] = forward_pass( ...
            x0, xtraj, utraj, ktraj, Ktraj, N, dt, param, Q, R, Qf, xd, ...
            J, DYNAMICS, xtraj0, utraj0);

        % Backward pass: compute local gains (k,K) around the NEW nominal
        [Ktraj, ktraj] = backward_pass( ...
            xtraj, utraj, ktraj, Ktraj, Q, R, Qf, xd, param, N, dt, ...
            GRADIENTS, xtraj0, utraj0);

        % Visualization of nominal trajectory evolution
        figure(1)
        plot(xtraj(1,:),xtraj(2,:))

        % Stop when cost improvement saturates
        if abs(J-Jlast) < 0.001
            break;
        end
        Jlast = J;
    end

    figure(1)
    plot(xtraj(1,:),xtraj(2,:))
    figure(2)
    plot(xtraj(5,:),xtraj(6,:))
end

%==========================================================================
% Q_terms: local Q-function derivatives for the backward pass
%
% Given stage cost derivatives (g*) and linearized dynamics (fx, fu),
% and value function derivatives (Vx, Vxx), form:
%   Qx, Qu, Qxx, Qux, Quu
%
% Also forms “regularized” (bar) versions Quxbar, Quubar by adding rho*I
% to Vxx inside the control Hessian terms (a Levenberg / trust-region knob).
%==========================================================================
function [Qx, Qu, Qxx, Qux, Quu, Quubar, Quxbar]= Q_terms ( ...
    gx, gu, gxx, gux, guu, fx, fu, Vx, Vxx)

    Qx  = gx  + fx'*Vx;
    Qu  = gu  + fu'*Vx;

    Qxx = gxx + fx'*Vxx*fx;
    Qux = gux + fu'*Vxx*fx;
    Quu = guu + fu'*Vxx*fu;

    % regularization on Vxx (rho) to stabilize Quu inversion
    rho = 1e-6;
    Quxbar = gux + fu'*(Vxx+rho*eye(length(Vxx)))*fx;
    Quubar = guu + fu'*(Vxx+rho*eye(length(Vxx)))*fu;
end

%==========================================================================
% gains: compute local policy δu = k + K δx
%==========================================================================
function [K, k]= gains(Qx,Qu,Qxx,Qux,Quu)
    Quu_reg = Quu + 1e-4*eye(size(Quu));
    k = -(Quu_reg \ Qu);
    K = -(Quu_reg \ Qux);
end

%==========================================================================
% Vterms: propagate value function derivatives backward
%==========================================================================
function [Vx, Vxx] = Vterms(Qx,Qu,Qxx,Qux,Quu,K,k)
    Vx  = Qx  + K'*Qu + Qux'*k + K'*Quu*k;
    Vxx = Qxx + K'*Qux + Qux'*K + K'*Quu*K;
    Vxx = 0.5*(Vxx + Vxx');  
end

%==========================================================================
% backward_pass
%
% Linearize discrete-time dynamics around nominal:
%   x_{k+1} ≈ x_k + f(x_k,u_k) dt
% so:
%   Fx ≈ I + (∂f/∂x) dt,   Fu ≈ (∂f/∂u) dt
%
% Then compute (k_k, K_k) backwards from terminal value.
% The stage cost here is a quadratic shaping about the previous nominal
% (xtraj00, utraj00): see cost_gradients().
%==========================================================================
function [Ktraj, ktraj] = backward_pass( ...
    xtraj, utraj, ktraj, Ktraj, Q, R, Qf, xd, param, N, dt, ...
    GRADIENTS, xtraj00, utraj00)

    % Terminal derivatives: ℓ_f(x_N)
    [gx, gu, gxx, gux, guu] = final_cost_gradients(xtraj(:,N),utraj(:,N-1)*0,xd,Qf);
    Vxx = gxx;
    Vx  = gx;

    for i=N-1:-1:1
        % Stage derivatives of shaping term about (x̄_i, ū_i)
        [gx, gu, gxx, gux, guu] = cost_gradients( ...
            xtraj(:,i),utraj(:,i),xd,Q,R,xtraj00(:,i),utraj00(:,i));

        % Continuous-time Jacobians at (x_i, u_i)
        [fx, fu]=GRADIENTS(0,xtraj(:,i),utraj(:,i),param);

        % Discretize: x_{k+1} = x_k + f dt
        fu = fu*dt;
        fx = (fx*dt+eye(length(fx)));

        % Q-function derivatives and local gains
        [Qx, Qu, Qxx, Qux, Quu,Quubar,Quxbar]= Q_terms ( ...
            gx, gu, gxx, gux, guu, fx, fu, Vx, Vxx);

        [Ktraj(:,:,i), ktraj(:,i)]= gains(Qx,Qu,Qxx,Quxbar,Quubar);

        % Value function update (uses unbarred Qux/Quu in this implementation)
        [Vx, Vxx] = Vterms (Qx,Qu,Qxx,Qux,Quu,Ktraj(:,:,i), ktraj(:,i));
    end
end

%==========================================================================
% cost: quadratic shaping term around (x0,u0) = (x̄_k, ū_k)
%==========================================================================
function J = cost(x,u,Q,R,x0,u0)
    J = 0.5*(x-x0)'*Q*(x-x0) + 0.5*(u-u0)'*R*(u-u0);
end

%==========================================================================
% final cost: terminal objective about xd
%==========================================================================
function Jf = final_cost(x,u,xd,Qf)
    Jf = 0.5*(x-xd)'*Qf*(x-xd);
end

%==========================================================================
% final cost gradients
%==========================================================================
function [gx, gu, gxx, gux, guu] = final_cost_gradients(x,u,xd,Qf)
    gx  = Qf*(x-xd);
    gu  = 0;
    gxx = Qf;
    gux = xd'*0;
    guu = 0;
end

%==========================================================================
% cost gradients for shaping term
%==========================================================================
function [gx, gu, gxx, gux, guu] = cost_gradients(x,u,xd,Q,R,x0,u0)
    gx  = Q*(x-x0);
    gu  = R*(u-u0);
    gxx = Q;
    gux = zeros(length(u),length(x));
    guu = R;
end

%==========================================================================
% forward_pass
%
% Given nominal (xtraj0, utraj0) and local gains (k,K), roll out:
%   u_i = u_i^{nom} + α k_i + K_i (x_i - x_i^{nom})
% with a simple backtracking line-search on α until cost decreases.
% After convergence, execution typically uses:
%   u_i = u_i^{nom} + K_i (x_i - x_i^{nom})
% (feedforward k is only used to update the nominal during planning).
%==========================================================================
function [xtraj, utraj, J] = forward_pass( ...
    x0, xtraj0, utraj0, ktraj, Ktraj, N, dt, param, Q, R, Qf, xd, ...
    J0, DYNAMICS, xtraj00, utraj00)

    J = 1e7;
    alpha = 10;

    % Backtracking until we achieve an improvement over J0
    while(J0<J)
        t = 0;
        x = x0;
        J = 0;

        for i=1:N-1
            xtraj(:,i) = x;

            % Forward-pass nominal update:
            %   u_i^{new} = u_i^{nom} + α k_i + K_i (x_i^{new} - x_i^{nom})
            u = utraj0(:,i) + alpha*ktraj(:,i) + Ktraj(:,:,i)*(x-xtraj0(:,i));
            utraj(:,i) = u;

            % Shaping cost around previous nominal (x̄_i, ū_i)
            J = J + cost(x,utraj(:,i),Q,R,xtraj00(:,i),utraj00(:,i));

            % Rollout dynamics (Euler)
            xdot = DYNAMICS(t,x,utraj(:,i),param);
            t = t + dt;
            x = x + xdot*dt;
        end

        xtraj(:,N) = x;

        % Terminal objective
        J = J + final_cost(x,u,xd,Qf);

        J0
        alpha = alpha/2
    end
end
