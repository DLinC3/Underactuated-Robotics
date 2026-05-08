% =========================================================================
% Cartpole Backpropagation Through Time Trajectory Optimization
% =========================================================================
function  cartpole_bptt_extend()

close all
clear all

% set up the dynamics
addpath("dynamical_systems/")
param=init_cartpole_params();
GRADIENTS=@cartpole_grads;
DYNAMICS=@cartpole_dynamics;
DRAW=@draw_cartpole;
rmpath("dynamical_systems/")

param.dt=1/20; %time step

param.x0=[0 0 0 0]'; % initial conditions
param.xd = [0 pi 0 0]'; % final desired state
param.iter=200; % number of iterations

param.Q=diag([5, 30, 1, 15]); % Cost on glider trajectory
param.R=0.01; % actuator cost
param.Qend=diag([2, 467.129702, 80, 20]); % final value cost

param.eta=0.005; % optimization update parameter
param.x0=[0 0 0 0]';

%bounds for plotting
param.xL=-2;
param.xU=2;
param.yL=-2*pi;
param.yU=2*pi;

dt=param.dt;%time step

N=50;%number of time steps
param.N=N;

utape0=0*ones(1,N);%initialize the initial open loop policy

ti=0;
t=ti:dt:(N-1)*dt;%time values

x0=param.x0;%set the initial condition

utape=bptt(utape0,x0,DYNAMICS,GRADIENTS,param);% run the bptt algorithm     
utape_extended = [utape, zeros(1, N)];  

%simulate the dynamics
xtape = zeros(4, 2*N);  
x = x0;
for k = 1:(2*N)
    xtape(:,k) = x;
    if k <= N
        u = utape(k);  
    else
        u = 0;         
    end
    xdot = DYNAMICS(t, x, u, param);
    x = x + xdot * dt;
end
%plot the resulting trajectory             
figure()
plot(xtape(1,:), xtape(2,:),0,pi,'o',0,0,'*')
axis([param.xL param.xU param.yL param.yU])

figure()
plot(utape_extended);

% draw the robot
for i = 1:size(xtape, 2)
    DRAW((i-1)*dt, xtape(:,i), param);
end

end