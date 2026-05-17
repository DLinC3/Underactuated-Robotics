
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function explores several sample-based trajectory optimization
% approaches for the perching glider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glider_gradient_free
    close all
    clear all

    addpath("dynamical_systems/")
    param=init_glider_params();
    DYNAMICS=@glider;
    DRAW=@draw_glider;
    rmpath("dynamical_systems/")
    
    %
    mode = 'ce';

    T = 0.75; %time horizon
    dt = 0.01;%time step
    N = floor(T/dt);% number of time steps
    nX = 7;%number of states
    nU = 1;%number of inputs
    
    t = zeros(1,N);%time
    x = zeros(nX,N);%state
    u = zeros(nU,N);%input

    % cost matrices
    Q = zeros(nX);
    Qf = eye(nX);
    Qf(1,1)=100;
    Qf(2,2)=100;
    Qf(3,3)=10;
    Qf(4,4)=10;
    Qf(5,5)=10;
    R = 0.1 * eye(nU);

    %desired final state
    xd = zeros(nX, 1);
    xd(1) = 0;
    xd(2) = 0;
    xd(3) = pi/4.0;
    xd(5) = 2;
    xd(6) = -2;
    param.xd = xd;
    
    %initial state
    x(:,1)=[-3.5;0.1;0;0;7;0;0];

    %switch between different gradient free methods
    switch mode
        case 'mppi'
            U=MPPI(Q,R,Qf,x(:,1),xd,N,dt,param, DYNAMICS);
        case 'ce'
            U=cross_entropy(Q,R,Qf,x(:,1),xd,N,dt,param, DYNAMICS);
    end
    
    %simulate the glider system
    for k =1:N-1
        x(:,k+1)=x(:,k)+ (DYNAMICS(t(k),x(:,k),U(:,k),param))*dt;
        t(k+1) = t(k)+dt;
    end
    
    %draw the glider system
    for k =1:N
        DRAW(t(k),x(:,k),param);
        drawnow;
    end  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The cross-entropy method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = cross_entropy(Q,R,Qf,x0,xd,N,dt,param, DYNAMICS)

num_samples = 30; %select the number of samples per iteration
num_elites = 10; %select the number of elite samples
iterations = 100; %select the number of total iterations
Jsave=[];
ksave=[];

U=zeros(length(R),N); %mean for the parameters
x = zeros(length(xd),N);
sigma2 = .1*ones(1,N); %choose the initial variance for parameter distribution
Jlast = 1e6;
for i=1:iterations
    i
    figure(7)
    cla;
    J = [];
    udu = zeros(size(U,1), N, num_samples); % store
    %loop through trajectory rollouts 
    for j=1:num_samples
        x(:,1) = x0;
        du = randn(size(U)) .* sqrt(sigma2); %sample from the parameter distribution
        dU(:,:,j) = du;
        u = U + du;
        udu(:,:,j) = u;
        t = 0;
        %compute the cost for the sampled parameters (FILL IN)
        [Jk,~,x] = sampleTrajectoryCosts(x0,xd,u,Q,R,Qf,dt,N,param,DYNAMICS);
        %stack the costs for each rollout
        J = [J;Jk];
        figure(7)
        hold on
        plot(x(1,:),x(2,:),'b')
        hold off        
    end
    
    %select the elite samples (FILL IN)
    % [Je, Ie]= ;  

    %find the sample with the minimum cost
    % [Jmin, Imin] = min(J);    
    
    % %select the parameters with minimum cost
    % Umin = udu(:,:,Imin); 
    % %select the elite parameters (FILL IN)
    % udue = ;

    [~, sorted_indices] = sort(J);% Sort by cost in ascending order
    elite = sorted_indices(1:num_elites);
    udu_elite = udu(:,:,elite);
    U_new = zeros(size(U));
    sigma2_new = zeros(size(sigma2));

    

    
    %use maximum likelihood esitmation to compute the new parameter
    %distribution (FILL IN)
    for j = 1:N
        % Calculate the elite sample mean and variance
        controls_at_k = udu_elite(:,j,:);   % Control input of all elite samples at time k
        U_new(:,j) = mean(controls_at_k, 3);% Update mean
        
        if num_elites > 1
            % Unbiased variance estimate (denominator is n-1)
            sigma2_new(j) = var(controls_at_k(:), 0); 
        else
            % to avoid zero variance
            sigma2_new(j) = sigma2(j); 
        end
    end
    sigma2_new = max(sigma2_new, 0.01);     % Set the minimum variance lower limit

    % rollout the mean for the new parameter distribution and evaluate the
    % cost
    U = U_new;
    sigma2 = sigma2_new;

    [~, Jk, x] = sampleTrajectoryCosts(x0,xd,U,Q,R,Qf,dt,N,param,DYNAMICS);
    Jc = sum(Jk);
    % Sig2 = diag(sigma2);
    % det(Sig2)   

    % U= Unew;
    % sigma2 = sigma2new;
    Jsave = [Jsave Jc];
    ksave = [ksave i];
end

figure(8)
plot(ksave,Jsave)
xlabel('iterations')
ylabel('Cost')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model predictive path-integral control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U=MPPI(Q,R,Qf,x0,xd,N,dt,param, DYNAMICS)

num_samples = 30; %select number of samples
iterations = 100; %select number of iterations
lambda = 10; % choose the lambda cost parameters
nu = 10;% set the variance of the uncertainty (FILL IN)
rho=10000; %set the scale factor

Jsave=[];
ksave=[];

J = [];
expS=[];
U=zeros(length(R),N);
u = zeros(length(R),N);
x = zeros(length(xd),N);

for i=1:iterations
    i
    figure(7)
    cla;
    J = [];
    expS=[];
    for j=1:num_samples
        x(:,1) = x0;
        %sample the input disturbance assuming a Gaussian noise (FILL IN)
        du = (sqrt(nu)/sqrt(rho)) * (randn(1, N)/sqrt(dt));
        dU{j} = du;
        u = U + du;
        t = 0;
        %roll out trajectories and evaluate the cost
        [~,Jk,x] = sampleTrajectoryCosts(x0,xd,u,Q,R,Qf,dt,N,param,DYNAMICS);
        figure(7)
        hold on
        plot(x(1,:),x(2,:),'b')
        hold off
        %compute the cost-to-gos for each time step
        Ji = getCostToGo(Jk);
        J=[J;Ji];
        expSi = exp(-Ji/lambda);%take the exponential of the cost-to-gos (FILL IN)
        expS=[expS;expSi];
    end
    for j=1:N
        ss= 0;
        su = 0;
        %compute the reward-weighted average MPPI iterative update law
        %(FILL IN)
        for s = 1:num_samples
            ss = ss + expS(s,j) * dU{s}(j);
            su = su + expS(s,j);
        end
        U(:,j) = U(:,j) + ss/su;
    end
    u = U;
    Jk = sampleTrajectoryCosts(x0,xd,u,Q,R,Qf,dt,N,param,DYNAMICS);
    Jc = sum(Jk);

    Jsave = [Jsave Jc];
    ksave = [ksave i];

end

figure(8)
plot(ksave,Jsave)
xlabel('iterations')
ylabel('Cost')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roll out the trajectories and sample from the costs (FILL IN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, Jk,x] = sampleTrajectoryCosts(x0,xd,u,Q,R,Qf,dt,N,param,DYNAMICS)
    x(:,1) = x0;
    t = 0;
    for k= 1:N
        Jk(k) = cost(x(:,k), u(:,k), Q, R, dt); %cost at each time step
        x(:,k+1) = x(:,k) + DYNAMICS(t, x(:,k), u(:,k), param) * dt; %integrate the one step dynamics
        t = t+dt;
    end
    Jk(N+1) = final_cost(x(:,N+1), xd, Qf); %compute the final cost
    J = sum(Jk);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute running cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = cost(x,u,Q,R,dt)

J = x'*Q*dt*x+u'*R*dt*u;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute final cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = final_cost(x,xd,Qf)
J = (x-xd)'*Qf*(x-xd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the cost-to-go for each time step from the instantaneous cost values Jk (FILL IN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = getCostToGo(Jk)
    S = zeros(size(Jk));
    S(end) = Jk(end);
    for k = length(Jk)-1:-1:1
        S(k) = Jk(k) + S(k+1);
    end
end
