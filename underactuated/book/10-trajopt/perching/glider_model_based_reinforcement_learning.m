%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   plan (iLQR) -> execute (sim) -> identify (ridge LS) -> re-plan
%
%   Learned planner dynamics use a linear-in-parameters aerodynamic residual:
%       xdot = f(x,u) + g(x,u) θ,
%       θ = [θ_Cl; θ_Cd; θ_Cm]
% where each coefficient is represented by RBF features φ(·) with bias:
%       C(·) = θ^T φ(·).
%
% Planner/simulator mismatch:
%   - Planner model uses features of (α, δ_e)
%   - Simulator uses features of (α, δ_e, V) (higher fidelity / different inputs)
%
% At iteration j:
%   1) iLQR plans on the current learned model (param.theta*), producing:
%        nominal (xtraj,utraj) and gains (Ktraj, ktraj)
%   2) Execute on simulator with uncertainty (random initial speed, data noise):
%        u_k = k_k + K_k (x_k - x_k^{nom})   (tracking law)
%   3) Identify θ via ridge LS from rollout data:
%        Φ := stack g(x_i,u_i) on affected components,
%        y := (xdot_data - f),
%        θ = (Φ'Φ + γI)^{-1} Φ' y
%   4) Write θ back into planner model and repeat.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glider_model_based_reinforcement_learning

    close all 
    clear all

    % Horizon / discretization
    T  = .75;
    dt = 0.01;
    N  = floor(T/dt);

    nX = 7;
    nU = 1;
    
    t = zeros(1,N);
    x = zeros(nX,N);
    u = zeros(nU,N);
    
    % iLQR needs dynamics gradients (∂f/∂x, ∂f/∂u)
    GRADIENTS=@glider_grads_7d_gaussian;

    % Planner vs simulator dynamics
    SIM_DYNAMICS=@glider_7d_gaussian;
    DYNAMICS=@glider_7d_gaussian;

    % Parameters and drawing
    param=init_glider_params();
    DRAW=@draw_glider;
 
    % Initial condition (nominal)
    x0= [-3.5;0.1;0;0;7;0;0];
    
    % iLQR nominal trajectories and gains
    xtraj=zeros(nX,N);
    utraj=zeros(nU,N);
    ktraj= zeros(nU,N);
    Ktraj = zeros(nU,nX,N);

    % Planning cost: terminal objective + shaping around previous nominal
    Q  = .1*eye(nX,nX);
    Qf = eye(nX,nX);
    Qf(1,1)=1000;
    Qf(2,2)=1000;
    Qf(3,3)=100;
    Qf(4,4)=1;
    Qf(5,5)=1;
    Qf(6,6)=1;
    R  = .1 * eye(nU);

    xd = [0;0;pi/4;0;2;-2;0];
    param.xd = xd;

    % Nominal wind setup
    param.use_wind = false;
    param.vx0= 0;
    param.vz0 =0;

    % ---------------------------------------------------------------------
    % Load learned RBF bases for planner model (μ, Σ) and initial θ
    % (Then we overwrite θ to zeros to start learning from scratch.)
    % ---------------------------------------------------------------------
    addpath("data/")
    coeffs=load('model_coeffs.mat');
    param.thetaCl=coeffs.thetaCl;
    param.thetaCd=coeffs.thetaCd;
    param.thetaCm=coeffs.thetaCm;
    param.muCl=coeffs.muCl;
    param.sigmaCl=coeffs.sigmaCl;
    param.muCm=coeffs.muCm;
    param.sigmaCm=coeffs.sigmaCm;
    param.muCd=coeffs.muCd;
    param.sigmaCd=coeffs.sigmaCd;
    param.mu=coeffs.muCl;
    param.sigma=coeffs.sigmaCl;

    sim_param = param;

    % ---------------------------------------------------------------------
    % Load higher-fidelity simulator coefficients (different feature inputs:
    % (α,δ_e,V) vs (α,δ_e) used by planner) -> deliberate mismatch.
    % ---------------------------------------------------------------------
    sim_coeffs=load('model_coeffs_vel.mat');

    sim_param.thetaCl=sim_coeffs.thetaCl;
    sim_param.thetaCd=sim_coeffs.thetaCd;
    sim_param.thetaCm=sim_coeffs.thetaCm;
    sim_param.muCl=sim_coeffs.muCl;
    sim_param.sigmaCl=sim_coeffs.sigmaCl;
    sim_param.muCm=sim_coeffs.muCm;
    sim_param.sigmaCm=sim_coeffs.sigmaCm;
    sim_param.muCd=sim_coeffs.muCd;
    sim_param.sigmaCd=sim_coeffs.sigmaCd;
    sim_param.mu=sim_coeffs.muCl;
    sim_param.sigma=sim_coeffs.sigmaCl;

    % Start planner with zero residual parameters (θ = 0)
    param.thetaCl=0*coeffs.thetaCl;
    param.thetaCd=0*coeffs.thetaCd;
    param.thetaCm=0*coeffs.thetaCm;

    Np = length(param.thetaCl);
    theta_hat = zeros(3*Np,1);

    % Learning loop
    num_iterations = 50;

    % Dataset accumulators
    xdata = [];
    udata = [];
    xdot_data = [];

    param.x0 = x0;
    param.plot = 1;
    
    rmpath("data/")
    for j=1:num_iterations

        % -------------------------------------------------------------
        % (1) plan (iLQR): optimize on current learned planner model
        % -------------------------------------------------------------
        [xtraj, utraj, ktraj, Ktraj] = iterative_LQR_reg( ...
            x0, xtraj, utraj, ktraj, Ktraj, N, dt, param, Q, R, Qf, xd, DYNAMICS, GRADIENTS);

        % -------------------------------------------------------------
        % (2) execute (sim): track nominal on simulator under uncertainty
        %     using standard tracking law (no α here):
        %         u_k = k_k + K_k (x_k - x_k^{nom})
        % -------------------------------------------------------------
        vmax = 7.5;
        vmin = 6.5;
        x(:,1)=[-3.5;0.1;0;0;vmin + (vmax-vmin).*rand(1,1);0;0];

        for k =1:N-1
            u(:,k) = utraj(:,k) + Ktraj(:,:,k)*(x(:,k)-xtraj(:,k));
            x(:,k+1)=x(:,k)+ (SIM_DYNAMICS(t(k),x(:,k),u(:,k),sim_param))*dt;
            t(k+1) = t(k)+dt;
        end
        
        % Performance metric: terminal objective J_j = ℓ_f(x_N)
        J = final_cost(x(:,N),u(:,N-1),xd,Qf);
        Jsave(j) = J;

        % Visualize execution
        for k =1:N
            DRAW(t(k),x(:,k),param);
            drawnow;
        end

        % Add measurement/process noise to data (simulation-only)
        x = x + 0.01*randn(nX,N);

        % -------------------------------------------------------------
        % (3) identify (ridge LS): estimate θ from rollout data
        %
        % Build finite-difference derivatives:
        %   xdot_data ≈ (x_{k+1}-x_k)/dt
        % and exploit linear structure:
        %   xdot = f(x,u) + g(x,u) θ
        % so y := xdot_data - f, Φ := g, solve ridge LS.
        % -------------------------------------------------------------
        theta_hat_last = theta_hat;

        xdata    = [xdata x(:,1:end-1)];
        udata    = [udata u(:,1:end-1)];
        xdot     = (x(:,2:end) - x(:,1:end-1))/dt;
        xdot_data=[xdot_data xdot];

        theta_hat = least_squares(t, xdata, udata, xdot_data, DYNAMICS, param, dt);

        % Write θ back into planner model
        param.thetaCl = theta_hat(1:Np);
        param.thetaCd = theta_hat(Np+1:2*Np);
        param.thetaCm = theta_hat(2*Np+1:3*Np);

        % -------------------------------------------------------------
        % (4) re-plan: next iteration uses updated planner parameters
        % -------------------------------------------------------------
        etheta = norm(theta_hat-theta_hat_last);
        etheta_save(j) = etheta;
        js(j) = j;

        figure(7)
        plot(js,etheta_save)
        xlabel('iteration')
        ylabel('norm of parameter difference')
        title('Parameter Change vs. Iteration')

        figure(8)
        plot(js, Jsave)
        xlabel('iteration')
        ylabel('final cost')
        title('Cost vs. Iteration')
    end
end


%==========================================================================
% least_squares
%
% Identification exploiting xdot = f(x,u) + g(x,u)θ:
%   - For each data point i:
%       compute (f,g) under param with θ temporarily set to zero
%       define y_i := xdot_data - f
%       stack g into Φ
%   - Solve ridge LS: θ = (Φ'Φ + γI)^{-1} Φ' y
%
% This implementation only uses aerodynamically affected components (5:7).
%==========================================================================
function theta_hat = least_squares(t, xdata, udata, xdot_data, DYNAMICS, param, dt)
    Phi=[];
    y = [];
    for i=1:length(xdata)
        x = xdata(:,i);
        u = udata(:,i);

        % Use measured/finite-difference derivative (only components 5:7)
        f0 = xdot_data(:,i);
        f0 = f0(5:7);

        % Evaluate nominal f and sensitivity g with θ=0 so that:
        %   xdot - f = g θ
        param.thetaCl=0*param.thetaCl;
        param.thetaCd=0*param.thetaCd;
        param.thetaCm=0*param.thetaCm;

        [f, g]= DYNAMICS(t,x,u,param);
        f = f(5:7);
        g = g(5:7,:);

        Phi = [Phi; g];
        y   = [y; f0-f];
    end

    gamma = 0.1;
    theta_hat = (Phi'*Phi + gamma*eye(size(Phi,2)))\(Phi'*y);
end



%==========================================================================
% cost
%==========================================================================
function J = cost(x,u,Q,R)
    J = 0.5*x'*Q*x + 0.5*u'*R*u;
end

%==========================================================================
% final cost: terminal objective
%==========================================================================
function Jf = final_cost(x,u,xd,Qf)
    Jf = 0.5*(x-xd)'*Qf*(x-xd);
end


%==========================================================================
% init_glider_params: baseline physical parameters for the learned model
%==========================================================================
function param=init_glider_params()

    % x,z,theta,phi_dot,x_dot,z_doth,theta_dot
    param.x0=[-3.4 0.1 0 0 7.4 0 0]';

    param.dt=1/119;
    param.g=9.81;
    param.rho= 1.292;

    dihedral = deg2rad(0);
    cm=0;

    param.S_w = 0.05*cos(dihedral) + 0.035 + 0.0035;
    param.S_e = 0.0147;
    param.l_w = cm;
    param.l_e = 0.022;
    param.l_h = 0.27-cm;
    param.m = 0.12;
    param.I=0.0015;

    param.cg_x=-0.065;
    param.cg_z=-0.025;
    param.vicon_offset=-0.052;
end


% =========================================================================
% glider_7d_gaussian
%
% Nonlinear glider model with linear-in-parameters aerodynamic residual:
%   - Computes RBF features φ(·) for Cl, Cd, Cm
%   - Forms coefficients: Cl = θ_Cl^T φ_Cl, etc.
%   - Returns:
%       xdot : state derivatives
%       F    : Jacobian w.r.t. θ (i.e., g(x,u) in xdot = f + gθ)
%
% Planner/simulator mismatch can be introduced by changing feature inputs:
%   - If sigma has length 3, features use [α; δ_e; V]
%   - Else features use [α; δ_e]
% =========================================================================
function [xdot, F] = glider_7d_gaussian(t,x,u,param)

m=param.m;
g=param.g;
rho=param.rho;
S_w=param.S_w;
S_e=param.S_e;
I=param.I;
l_h=param.l_h;
l_w=param.l_w;
l_e=param.l_e;

Sw=S_w;
Se=S_e;
l=l_h;
lw=l_w;
le=l_e;

muCl=param.muCl; sigmaCl=param.sigmaCl;
muCd=param.muCd; sigmaCd=param.sigmaCd;
muCm=param.muCm; sigmaCm=param.sigmaCm;

thetaCl=param.thetaCl;
thetaCd=param.thetaCd;
thetaCm=param.thetaCm;

q3 = x(3,:); q4 = x(4,:);
qdot10 = x(5,:); qdot20 = x(6,:); qdot3 = x(7,:); qdot4=0*u(1);

% Inject wind into inertial velocity components
qdot1=qdot10+param.vx0;
qdot2=qdot20+param.vz0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flat-plate-like forces (Fw, Fe) from wing/elevator geometry 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xwdot = qdot1 - lw*qdot3.*sin(q3);
zwdot = qdot2 + lw*qdot3.*cos(q3);
alpha_w = q3 - atan2(zwdot,xwdot);
Fw = rho*Sw*sin(alpha_w).*(zwdot.^2+xwdot.^2);

xedot = qdot1 + l*qdot3.*sin(q3) + le*(qdot3+qdot4).*sin(q3+q4);
zedot = qdot2 - l*qdot3.*cos(q3) - le*(qdot3+qdot4).*cos(q3+q4);
alpha_e = q3+q4-atan2(zedot,xedot);
Fe = rho*Se*sin(alpha_e).*(zedot.^2+xedot.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learned aerodynamic augmentation via RBF features (linear in θ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = q3 - atan(qdot2/qdot1);      % AoA proxy
V_norm = sqrt(qdot1^2+qdot2^2);      % speed norm

n_L1=-qdot2/V_norm; n_L2= qdot1/V_norm;  % lift direction
n_D1=-qdot1/V_norm; n_D2=-qdot2/V_norm;  % drag direction

q=0.5*rho*V_norm^2*S_w;              % dynamic pressure (scaled)

% Feature inputs: either [α; δ_e] or [α; δ_e; V]
if length(sigmaCl) == 3
    phiCl= gaussian([alpha q4 V_norm]',muCl,sigmaCl);
    phiCd= gaussian([alpha q4 V_norm]',muCd,sigmaCd);
    phiCm= gaussian([alpha q4 V_norm]',muCm,sigmaCm);
else
    phiCl= gaussian([alpha q4]',muCl,sigmaCl);
    phiCd= gaussian([alpha q4]',muCd,sigmaCd);
    phiCm= gaussian([alpha q4]',muCm,sigmaCm);
end

% Coefficients are linear in θ
Cl=thetaCl'*phiCl;
Cd=thetaCd'*phiCd;
Cm=thetaCm'*phiCm;

% Lift/drag/moment contributions
L = q*Cl*[n_L1;n_L2];
D = q*Cd*[n_D1;n_D2];
M = q*Cm;

% g(x,u): sensitivity of xdot(5:7) w.r.t θ (stacked for LS identification)
F = zeros(length(x),length(thetaCl)+length(thetaCd)+length(thetaCm));
F(5:7,:) = [n_L1*q*phiCl'/m n_D1*q*phiCd'/m phiCm'*0;
            n_L2*q*phiCl'/m n_D2*q*phiCd'/m phiCm'*0;
            phiCl'*0       phiCd'*0       q*phiCm'/I];

F1=L+D;

% Assemble xdot
xdot(4,:)=u(1);
xdot(5:6,:) = (sum(F1,2))/m;
xdot(7,:)=M/I;

% Add baseline flat-plate forces
xdot(5,:)=xdot(5,:)+(-(Fw.*sin(q3) + Fe.*sin(q3+q4))./m);
xdot(6,:)=xdot(6,:)+((Fw.*cos(q3) + Fe.*cos(q3+q4)-m*g)./m);
xdot(7,:)=xdot(7,:)+((Fw.*lw - Fe.*(l*cos(q4)+le))./I);

% Kinematics
xdot(1,:)=x(5,:);
xdot(2,:)=x(6,:);
xdot(3,:)=x(7,:);
end

% ============================================
% gaussian: RBF feature map φ with bias
%
% φ_k(x) = exp(-1/2 (x-μ_k)^T Σ^{-1} (x-μ_k)),  plus a constant 1 feature.
% ============================================
function phi=gaussian(x,mu,sigma)
phi=[];
[M N]=size(mu);
for k=1:N
    p= exp(-0.5*(x-mu(:,k))'*inv(sigma)*(x-mu(:,k)));
    phi=[phi p];
end
phi=[phi 1];
phi=phi';
end


% ============================================
% draw_glider: visualization
% ============================================
function draw_glider(t,x,param)
xd=param.xd;

sc = 2;
lw = -0.15*sc;
lh = 0.45*sc;
le = 0.04*sc;
mac = .1145*sc;

ct = cos(x(3)); st = sin(x(3));
ctp = cos(x(3)+x(4)); stp = sin(x(3)+x(4));

figure(25); clf; hold on;

plot(x(1),x(2),'r.','MarkerSize',10);

if (lw < 0)
    line([x(1) x(1)-lh*ct],[x(2) x(2)-lh*st],...
        'LineWidth',1,'Color',[0 0 0]);
else
    line([x(1)+(lw+mac/2)*ct x(1)-lh*ct],...
        [x(2)+(lw+mac/2)*st x(2)-lh*st],...
        'LineWidth',1,'Color',[0 0 0]);
end

line([x(1)+(lw+mac/2)*ct x(1)+(lw-mac/2)*ct],...
    [x(2)+(lw+mac/2)*st x(2)+(lw-mac/2)*st],...
    'LineWidth',3,'Color',[0 .2 1]);

line([x(1)-lh*ct x(1)-lh*ct-2*le*ctp],...
    [x(2)-lh*st x(2)-lh*st-2*le*stp],...
    'LineWidth',3,'Color',[0 .2 1]);

axis equal; axis([-4 1 -1 1.75]);
title(['time: ',num2str(t,2),' s']);
xlabel('x (m)'); ylabel('z (m)');

plot(xd(1),xd(2),'ko','MarkerSize',5,'MarkerFaceColor',[0 0 0]);
end
