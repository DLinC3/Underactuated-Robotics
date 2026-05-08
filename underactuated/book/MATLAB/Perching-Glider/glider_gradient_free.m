%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-Order trajectory optimization for a perching glider (CEM / MPPI)
%
% Notes to My Future Self:
% - Both methods are "rollout dominated": repeatedly sample open-loop controls,
%   simulate forward with Euler integration, evaluate cost, then update U.
% - CEM: fit a diagonal Gaussian over u_k via elite-set MLE (mean/variance).
% - MPPI: treat samples as stochastic control perturbations and do a
%   reward-weighted (exp(-cost-to-go/lambda)) average of perturbations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glider_gradient_free
    close all
    clear all

    addpath("dynamical_systems/")
    param=init_glider_params();
    DYNAMICS=@glider;
    DRAW=@draw_glider;
    rmpath("dynamical_systems/")
    
    % Switch between 'ce' or 'mppi'
    mode = 'mppi'; % 

    % Horizon + discretization
    T = 0.75;      % time horizon
    dt = 0.01;     % time step
    N = floor(T/dt); % number of time steps
    nX = 7;        % number of states
    nU = 1;        % number of inputs
    
    % Preallocate for nominal rollout (for visualization only)
    t = zeros(1,N);
    x = zeros(nX,N);

    % Cost matrices (running + terminal)
    Q = zeros(nX);
    Qf = eye(nX);
    Qf(1,1)=100;
    Qf(2,2)=100;
    Qf(3,3)=10;
    Qf(4,4)=10;
    Qf(5,5)=10;
    R = 0.1 * eye(nU);

    % Desired final state
    xd = zeros(nX, 1);
    xd(1) = 0;
    xd(2) = 0;
    xd(3) = pi/4.0;
    xd(5) = 2;
    xd(6) = -2;
    param.xd = xd;
    
    % Initial state
    x(:,1)=[-3.5;0.1;0;0;7;0;0];

    % Optimize an open-loop sequence U = [u_1,...,u_N]
    switch mode
        case 'mppi'
            U=MPPI(Q,R,Qf,x(:,1),xd,N,dt,param, DYNAMICS);
        case 'ce'
            U=cross_entropy(Q,R,Qf,x(:,1),xd,N,dt,param, DYNAMICS);
    end
    
    % Simulate the glider system using the optimized open-loop sequence U
    for k =1:N-1
        x(:,k+1)=x(:,k) + DYNAMICS(t(k),x(:,k),U(:,k),param)*dt; % Euler step
        t(k+1) = t(k)+dt;
    end
    
    % Draw the glider system (no video recording)
    figure(1);
    for k =1:N
        DRAW(t(k),x(:,k),param);
        drawnow;
    end  
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Cross-Entropy Method (CEM)
%
% Maintains per-time-step Gaussian: u_k ~ N(U_k, sigma_k^2) (diagonal)
% Each iteration:
%  1) sample M trajectories
%  2) select elite set E (lowest cost)
%  3) update mean/variance by MLE on elites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = cross_entropy(Q,R,Qf,x0,xd,N,dt,param, DYNAMICS)

num_samples = 30;   % number of rollouts per iteration
num_elites  = 10;   % elite set size
iterations  = 100;  % number of CEM iterations

Jsave=[];
ksave=[];

U      = zeros(length(R),N); % current mean (open-loop sequence)
x      = zeros(length(xd),N);
sigma2 = .1*ones(1,N);       % per-time-step variance (diagonal)

for i=1:iterations
    disp(i);

    % Visualization of sampled rollouts
    figure(7); cla;
    title(sprintf('CEM samples, iteration %d', i));

    J = [];
    udu = zeros(size(U,1), N, num_samples); % store sampled controls u = U + du

    % Rollout loop: sample u-sequences and evaluate cost
    for j=1:num_samples
        x(:,1) = x0;

        % Sample perturbation du_k ~ N(0, sigma_k^2)
        du = randn(size(U)) .* sqrt(sigma2);

        u = U + du;               % sampled open-loop control sequence
        udu(:,:,j) = u;           % store for elite stats

        % Evaluate trajectory cost under this sampled sequence
        [Jk_scalar,~,x] = sampleTrajectoryCosts(x0,xd,u,Q,R,Qf,dt,N,param,DYNAMICS);

        % Stack costs for elite selection
        J = [J; Jk_scalar];

        % Plot rollout in (x,z) or (x(1), x(2)) plane
        hold on
        plot(x(1,:),x(2,:),'b')
        hold off
        drawnow;
    end
    
    % Elite selection: lowest-cost samples
    [~, sorted_indices] = sort(J);           % ascending
    elite = sorted_indices(1:num_elites);
    udu_elite = udu(:,:,elite);

    % MLE update of diagonal Gaussian parameters (mean + variance)
    U_new      = zeros(size(U));
    sigma2_new = zeros(size(sigma2));
    for k = 1:N
        controls_at_k = udu_elite(:,k,:);          % all elites at time k
        U_new(:,k) = mean(controls_at_k, 3);       % mean over elites

        if num_elites > 1
            % var(...,0) uses unbiased normalization (n-1)
            sigma2_new(k) = var(controls_at_k(:), 0);
        else
            sigma2_new(k) = sigma2(k);             % avoid zero variance
        end
    end

    % Prevent variance collapse (keeps exploration alive)
    sigma2_new = max(sigma2_new, 0.01);

    % Update distribution
    U      = U_new;
    sigma2 = sigma2_new;

    % Track cost of current mean control sequence
    [~, Jk_vec, ~] = sampleTrajectoryCosts(x0,xd,U,Q,R,Qf,dt,N,param,DYNAMICS);
    Jc = sum(Jk_vec);

    Jsave = [Jsave Jc];
    ksave = [ksave i];
    data_ce = [ksave(:), Jsave(:)];
    writematrix(data_ce, 'ce_cost.csv');
end

figure(8);
plot(ksave,Jsave);
xlabel('iterations');
ylabel('Cost');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Predictive Path-Integral control (MPPI)
%
% Uses stochastic control perturbations:
%   u_k^(j) = U_k + delta u_k^(j)
% with delta u_k sampled i.i.d. Gaussian each step.
%
% For each rollout j, compute cost-to-go S_k^(j), form weights
%   w_k^(j) ∝ exp( -S_k^(j) / lambda )
% and update
%   U_k ← U_k + (Σ_j w_k^(j) delta u_k^(j)) / (Σ_j w_k^(j)).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U=MPPI(Q,R,Qf,x0,xd,N,dt,param, DYNAMICS)

num_samples = 30;
iterations  = 100;

lambda = 10;   % temperature in exp(-S/lambda) (your text called it kappa)
nu     = 10;   % noise variance parameter
rho    = 10000;% noise scale (larger rho => smaller perturbations)

Jsave=[];
ksave=[];

U = zeros(length(R),N);     % nominal open-loop sequence
x = zeros(length(xd),N);

for i=1:iterations
    disp(i);

    figure(7); cla;
    title(sprintf('MPPI samples, iteration %d', i));

    J    = [];   % cost-to-go for each rollout (rows) and time (cols)
    expS = [];   % exp(-S/lambda) weights

    dU = cell(1,num_samples); % store sampled perturbations

    for j=1:num_samples
        x(:,1) = x0;

        % Sample white-noise-like perturbation with dt scaling:
        % du_k ~ N(0, nu/(rho*dt))  (matches your write-up)
        du = (sqrt(nu)/sqrt(rho)) * (randn(1, N)/sqrt(dt));
        dU{j} = du;

        u = U + du;

        % Rollout and evaluate instantaneous costs Jk (and final cost)
        [~,Jk,x] = sampleTrajectoryCosts(x0,xd,u,Q,R,Qf,dt,N,param,DYNAMICS);

        hold on
        plot(x(1,:),x(2,:),'b')
        hold off
        drawnow;

        % Convert instantaneous costs to cost-to-go S_k = Σ_{ℓ=k}^{N+1} J_ℓ
        Ji = getCostToGo(Jk);     % 1 x (N+1)
        J  = [J; Ji];

        % Path-integral weights (per time step)
        expSi = exp(-Ji/lambda);
        expS  = [expS; expSi];
    end

    % MPPI update: reward-weighted average of perturbations, per time step
    for k=1:N
        ss = 0;  % weighted sum of perturbations
        su = 0;  % sum of weights
        for s = 1:num_samples
            ss = ss + expS(s,k) * dU{s}(k);
            su = su + expS(s,k);
        end
        U(:,k) = U(:,k) + ss/su;
    end

    % Track cost of current mean control sequence (vector Jk, then sum)
    [~,Jk_vec,~] = sampleTrajectoryCosts(x0,xd,U,Q,R,Qf,dt,N,param,DYNAMICS);
    Jc = sum(Jk_vec);

    Jsave = [Jsave Jc];
    ksave = [ksave i];
    data_mppi = [ksave(:), Jsave(:)];
    writematrix(data_mppi, 'mppi_cost.csv');
end

figure(8);
plot(ksave,Jsave);
xlabel('iterations');
ylabel('Cost');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roll out trajectory under a given open-loop control sequence u_{1:N}
% and compute:
%  - J  : total cost (scalar)
%  - Jk : instantaneous costs + terminal cost (length N+1)
%  - x  : state trajectory (size nX x (N+1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, Jk,x] = sampleTrajectoryCosts(x0,xd,u,Q,R,Qf,dt,N,param,DYNAMICS)
    x(:,1) = x0;
    t = 0;

    for k= 1:N
        Jk(k) = cost(x(:,k), u(:,k), Q, R, dt); % running cost at step k
        x(:,k+1) = x(:,k) + DYNAMICS(t, x(:,k), u(:,k), param) * dt; % Euler
        t = t+dt;
    end

    Jk(N+1) = final_cost(x(:,N+1), xd, Qf); % terminal cost
    J = sum(Jk);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running cost: x'Qx*dt + u'Ru*dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = cost(x,u,Q,R,dt)
    J = x'*Q*dt*x + u'*R*dt*u;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Terminal cost: (x-xd)' Qf (x-xd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = final_cost(x,xd,Qf)
    J = (x-xd)'*Qf*(x-xd);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert instantaneous costs Jk into cost-to-go:
%   S(k) = Jk(k) + Jk(k+1) + ... + Jk(end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = getCostToGo(Jk)
    S = zeros(size(Jk));
    S(end) = Jk(end);
    for k = length(Jk)-1:-1:1
        S(k) = Jk(k) + S(k+1);
    end
end
