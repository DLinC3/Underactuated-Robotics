%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online identification in EKF:
%
% Nonlinear dynamics (7D glider) depend on:
%   - aerodynamic parameters θ = [l_w; S_e]
%   - wind disturbance w = [v_x^0; v_z^0] (treated as time-varying in sim)
%
% The loop is:
%   1) Design an LQR controller around a desired operating point (x_d, u_d)
%      using numerical linearization (A,B) ≈ (∂f/∂x, ∂f/∂u).
%   2) Fly closed-loop with the true nonlinear plant and time-varying wind.
%   3) Run an augmented-state EKF on \bar x = [x; θ] with constant-parameter
%      process model \dot θ = 0, and measurement y = H \bar x + v (H=[I 0]).
%
% Outputs:
%   - perching trajectory
%   - θ_hat(t) traces from EKF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glider_nonlinear_parameter_estimation
    close all
    clear all

    % Plant and estimator models:
    %   DYNAMICS:       true 7D plant dynamics f(x,u;param,w)
    %   DYNAMICS_EKF:   augmented dynamics \bar f([x;θ],u) = [f(x,u;θ); 0]
    param=init_glider_params();
    DYNAMICS=@pglider_7d;
    DYNAMICS_EKF=@pglider_7d_augmented;
    DRAW=@draw_glider;

    % Horizon
    T = 4;
    dt = 0.01;
    N = floor(T/dt);
    nX = 7;
    nU = 2;
    
    % Allocate
    t = zeros(1,N);
    x = zeros(nX,N);

    % ---------------------------------------------------------------------
    % LQR design around desired operating point (xd, ud≈0)
    % (Reduced state ignores x-position to focus on z/theta/... regulation.)
    % ---------------------------------------------------------------------
    Q = eye(nX);
    Q(2,2)=1000;           % emphasize z regulation
    R = 0.1 * eye(nU);

    xd = zeros(nX, 1);
    xd(1) = 0;
    xd(2) = 0;
    xd(3) = 0.05;
    xd(5) = 9;
    xd(6) = 0;
    param.xd = xd;

    % Nominal wind (will be varied during the simulation)
    param.vx0 = -1;
    param.vz0 = 0.5;

    % Numerical linearization at (xd, u=0) to build LQR gain
    u0 = zeros(nU,1);
    [A, B] = NumGrads(DYNAMICS,xd,u0,param);

    % LQR on reduced state (ignore x-position: drop state 1)
    [K, S] = lqr(A(2:end,2:end), B(2:end,:), Q(2:end,2:end), R);

    % Reduced desired state for the reduced feedback
    xd = xd(2:end);
    
    % ---------------------------------------------------------------------
    % Wind process (in sim): w_{k+1} = w_k + W w_k + noise
    % (EKF here does NOT estimate wind; it treats wind as part of param.)
    % ---------------------------------------------------------------------
    W = -0.001*eye(2);
    w(:,1)=[param.vx0;param.vz0];

    % Initial true state
    x(:,1)=[-3.5;0.5;0;0;7;0;0];

    % EKF augmented state: xhat = [x; l_w; S_e]
    xhat(:,1) = [x(:,1);0;0];

    % EKF covariance + noise models
    Pk  = eye(length(xhat(:,1))) * 0.1;                % initial covariance
    Rest =  0.01 * eye(nX);                            % measurement noise
    Qest = diag([0.01*ones(1,7), 0.001, 0.001]);       % process noise

    % ---------------------------------------------------------------------
    % Closed-loop simulation with online augmented-state EKF
    % ---------------------------------------------------------------------
    for k =1:N-1
        % LQR control law on reduced state (excluding x-position)
        U(:,k) = -K*(x(2:end,k) - xd);

        % Update the plant wind parameters for this step
        param.vx0 = w(1,k);
        param.vz0 = w(2,k);

        % True plant rollout (Euler)
        x(:,k+1)=x(:,k)+ (DYNAMICS(t(k),x(:,k),U(:,k),param))*dt;

        % EKF predict/update on augmented model
        [xhat(:,k+1), Pk] = extended_kalman_filter( ...
            t, x(:,k), xhat(:,k), U(:,k), Pk, ...
            DYNAMICS_EKF, @NumGrads, param, Qest, Rest);

        % Log parameter estimates θ_hat = [l_w; S_e]
        theta_hat(:,k) = xhat(end-1:end);

        % Wind evolves (simulation-only)
        w(:,k+1) = w(:,k) + W*w(:,k) + 0.01*randn(2,1);

        t(k+1) = t(k)+dt;
    end
    
    % Draw flight (no video recording)
    for k =1:N
        % draw_glider expects [x; wind] appended for visualization
        DRAW(t(k),[x(:,k);w(1,k);w(2,k)],param);
        drawnow;
    end  

    figure(1)
    plot(x(1,:),x(2,:))
    xlabel('x (m)')
    ylabel('y (m)')

    figure(2)
    plot(t(1:end-1),theta_hat')
    xlabel('t (s)')
    ylabel('parameters ')
end

%==========================================================================
% NumGrads: numerical Jacobians (∂f/∂x, ∂f/∂u) around (x0,u0)
%==========================================================================
function [dFdx, dFdu] = NumGrads(DYNAMICS,x0,u0,param)

N=length(x0);
M=length(u0);

dFdx=zeros(N,N);
dFdu=zeros(N,M);
delta=1e-5;

for i=1:N
    
    xU=x0;
    xL=x0;
    
    xU(i)=x0(i)+delta;
    xL(i)=x0(i)-delta;
        
    xdotU=DYNAMICS(0,xU,u0,param);
    xdotL=DYNAMICS(0,xL,u0,param);
    
    dFdx(:,i)=(xdotU-xdotL)/(xU(i)-xL(i));
end

for i=1:M
    
    uU=u0;
    uL=u0;
    
    uU(i)=u0(i)+delta;
    uL(i)=u0(i)-delta;
        
    xdotU=DYNAMICS(0,x0,uU,param);
    xdotL=DYNAMICS(0,x0,uL,param);
    
    dFdu(:,i)=(xdotU-xdotL)/(uU(i)-uL(i));
end
end

%==========================================================================
% extended_kalman_filter
%
% Augmented-state EKF for \bar x:
%   predict:  \bar x_{k+1|k} = \bar x_k + \bar f(\bar x_k, u_k) dt
%   linearize: A = ∂\bar f/∂\bar x
%   update with measurement y = H \bar x + v  (here H=[I 0])
%==========================================================================
function [xhat, Pk] = extended_kalman_filter(t, x, xhat, u, Pk, DYNAMICS, GRADIENTS, param, Qs, Rs)
    dt = param.dt; % time-step used inside the EKF model
    nY = length(param.xd); % number of measurements (here: 7 states)

    % Predict (Euler)
    xhat = xhat + DYNAMICS(t, xhat, u, param)*dt;

    I = eye(length(xhat));

    % Linearize the process model
    [A, B]= GRADIENTS(DYNAMICS, xhat, u, param);
    Ak = I + A*dt;

    % Covariance predict
    Pk = Ak*Pk*Ak' + Qs;

    % Measurement model: y = H xhat + v, H=[I 0]
    Hk = [eye(nY), zeros(nY, length(xhat) - nY)];

    % Innovation
    ek = x - Hk*xhat;
    Sk =  Hk*Pk*Hk' + Rs;

    % Kalman gain
    Kk = Pk*Hk'/Sk;

    % Update
    xhat = xhat + Kk*ek;
    Pk = (I - Kk*Hk)*Pk;
end

%==========================================================================
% init_glider_params: physical + geometric parameters for the glider
%==========================================================================
function param=init_glider_params()

    % x,z,theta,phi_dot,x_dot,z_doth,theta_dot
    param.x0=[-3.4 0.1 0 0 7.4 0 0]';
    
    param.dt=1/119; % time step (note: EKF uses this dt)
    param.g=9.81;   % m/s^2
    param.rho= 1.292; % kg/m^3
    
    dihedral = deg2rad(0); % dihedral angle
    cm=0;

    % Areas / lengths
    param.S_w = 0.05*cos(dihedral) + 0.035 + 0.0035; % wing + fus + tail
    param.S_e = 0.0147;
    param.l_w = -0.03+cm;
    param.l_e = 0.022;
    param.l_h = 0.27-cm;

    % Mass / inertia
    param.m = 0.12;
    param.I=0.0015;
     
    % geometry offsets
    param.cg_x=-0.065;
    param.cg_z=-0.025;
    param.vicon_offset=-0.052;
end

%==========================================================================
% pglider_7d_augmented
%
% Augmented process model for EKF:
%   \bar x = [x; θ],   θ = [l_w; S_e]
%   \dot{\bar x} = [ f(x,u;θ,w); 0 ]
%==========================================================================
function xdota = pglider_7d_augmented(t,xa,u,param)

    % Inject parameter estimates from augmented state into param struct
    param.l_w = xa(8);
    param.S_e = xa(9);
    
    % Extract physical state
    x = xa(1:7);
    
    % Nonlinear glider dynamics
    xdot = pglider_7d(t,x,u,param);
    
    % Constant-parameter model: \dot θ = 0
    xdota = zeros(9, 1);
    xdota(1:7) = xdot;
    xdota(8:9) = 0; 
end
    
%==========================================================================
% pglider_7d: 7D planar glider dynamics with wind
% State:
%   x = [x; z; theta; phi; xdot; zdot; thetadot]
% Control:
%   u(1) = phidot command
%   u(2) = thrust-like term (appears in xdot/zdot equations)
% Wind:
%   param.vx0, param.vz0 are injected each step (simulation loop)
%==========================================================================
function xdot = pglider_7d(t,x,u,param)
    Sw=param.S_w;
    Se=param.S_e;
    l=param.l_h;
    lw=param.l_w;
    le=param.l_e;
    
    m=param.m;
    g=param.g;
    rho=param.rho;
    I=param.I;
    
    vx0=param.vx0;
    vz0=param.vz0;
    
    bw=0; % (unused here, kept from original code)
    
    q1 = x(1,:); q2 = x(2,:);  q3 = x(3,:); q4 = x(4,:);
    qdot1 = x(5,:); qdot2 = x(6,:); qdot3 = x(7,:); qdot4=u(1);
    q5=0;q6=0;
     
    % Wing relative velocity -> wing AoA -> wing force
    xwdot = qdot1+bw*u(2)*cos(q3) - vx0 - q5 - lw*qdot3.*sin(q3);
    zwdot = qdot2+bw*u(2)*sin(q3) - vz0 - q6 + lw*qdot3.*cos(q3);
    alpha_w = q3 - atan2(zwdot,xwdot);
    Fw = rho*Sw*sin(alpha_w).*(zwdot.^2+xwdot.^2);
    
    % Elevator relative velocity -> elevator AoA -> elevator force
    xedot = qdot1+bw*u(2)*cos(q3) - vx0 - q5 + l*qdot3.*sin(q3) + le*(qdot3+qdot4).*sin(q3+q4); 
    zedot = qdot2+bw*u(2)*sin(q3) - vz0 - q6 - l*qdot3.*cos(q3) - le*(qdot3+qdot4).*cos(q3+q4);
    alpha_e = q3+q4-atan2(zedot,xedot);
    Fe = rho*Se*sin(alpha_e).*(zedot.^2+xedot.^2);
    
    xdot=x;
    xdot(1,:)=x(5,:);
    xdot(2,:)=x(6,:);
    xdot(3,:)=x(7,:);
    xdot(4,:)=u(1);
    
    xdot(5,:)=(-(Fw.*sin(q3) + Fe.*sin(q3+q4))./m)+u(2)*cos(q3)./m;
    xdot(6,:)=((Fw.*cos(q3) + Fe.*cos(q3+q4))./m-g)+u(2)*sin(q3)./m;
    xdot(7,:)=((Fw.*lw - Fe.*(l*cos(q4)+le))./I);
end

%==========================================================================
% draw_glider: visualization (uses appended wind components for quiver)
%==========================================================================
function draw_glider(t,x,param)
xd=param.xd;

% x: [x z theta phi dot(x,z,theta,phi) windx windz] for visualization
sc = 2;
lw = -0.15*sc;
lh = 0.45*sc;
le = 0.04*sc;
mac = .1145*sc;

ct = cos(x(3)); st = sin(x(3));
ctp = cos(x(3)+x(4)); stp = sin(x(3)+x(4));

figure(25); clf; hold on;

% Draw wind velocity field (quiver)
Xl=-4+x(1);
Xu=1+x(1);
Yu=1.0;
Yl=-1.0;
dx=0.5;
dy=0.5;
X=Xl:dx:Xu;
Y=Yl:dy:Yu;
[XX,YY] = meshgrid(X,Y);
[m n]=size(XX);

Us=x(8)*ones(m,n)+param.vx0;
Vs=x(9)*ones(m,n);
quiver(XX,YY,Us,Vs,'Color','blue');

% Center of gravity
plot(x(1),x(2),'r.','MarkerSize',10);

% Fuselage
if (lw < 0)
    line([x(1) x(1)-lh*ct],[x(2) x(2)-lh*st],...
        'LineWidth',1,'Color',[0 0 0]);
else
    line([x(1)+(lw+mac/2)*ct x(1)-lh*ct],...
        [x(2)+(lw+mac/2)*st x(2)-lh*st],...
        'LineWidth',1,'Color',[0 0 0]);
end

% Wing
line([x(1)+(lw+mac/2)*ct x(1)+(lw-mac/2)*ct],...
    [x(2)+(lw+mac/2)*st x(2)+(lw-mac/2)*st],...
    'LineWidth',3,'Color',[0 .2 1]);

% Elevator
line([x(1)-lh*ct x(1)-lh*ct-2*le*ctp],...
    [x(2)-lh*st x(2)-lh*st-2*le*stp],...
    'LineWidth',3,'Color',[0 .2 1]);

axis equal; axis([-4+x(1) 1+x(1) -1 1.75]);
title(['time: ',num2str(t,2),' s']);
xlabel('x (m)'); ylabel('z (m)');

% Perch marker
plot(xd(1),xd(2),'ko','MarkerSize',5,'MarkerFaceColor',[0 0 0]);
end
