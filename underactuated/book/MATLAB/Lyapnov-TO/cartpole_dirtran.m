%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function computes optimal swing-up trajectory via direct transcription
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cartpole_dirtran

clear all
close all

%set up the dynamics
addpath("dynamical_systems/")
param=init_cartpole_params();
DYNAMICS=@cartpole_dynamics;
DRAW=@draw_cartpole;
addpath("dynamical_systems/")


N = 100;% FIL IN the number of timesteps
param.dt= 0.05; %FILL IN the discrete timestep

%set up the optimization parameters
param=opt_param(param);
param.nX=4; %state size
param.nU=1; %input size
param.xd=[0 pi]; %desried final state

param.plot=0;
tic
% generate a trajeectory via direct transcription
[w,F,info] = dirtran2(N,DYNAMICS,param); 
toc

info

% convert optimization parameters into state and control trajectories
xtape = reshape(w(2:(1+param.nX*N)),param.nX,N);
utape = reshape(w(param.nX*N+2:end),param.nU,N);

%plot the trajectory in x1-x2 space
figure(2)
hold on
plot(xtape(1,:),xtape(2,:));
plot(0,0,'ro')
hold off

% plot the x1 state vs time
figure(3)
hold on
plot(xtape(1,:));
plot(0,0,'ro')
hold off

%plot the x2 state vs time
figure(4)
hold on
plot(xtape(2,:));
plot(0,0,'ro')
hold off

%plot the input trajectory vs time
figure(5)
plot(utape)

% draw the cartpole
figure(25)
for i=1:size(xtape,2)
    DRAW((i-1)*param.dt,xtape(:,i),param);
end  
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method sets up the optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function param=opt_param(param)

param.opt.x0=[0 0 0 0]'; %initial state
param.opt.xd=[0 pi 0 0]'; %final state
%FILL IN bounds on final state
err = 0.05;
param.opt.xfmax=[err; pi+0.2*err; err; err];
param.opt.xfmin=[-err; pi-0.2*err; -err; -err];

%bounds on state
param.opt.xmax=[10 10 10 10];
param.opt.xmin=[-10 -10 -10 -10];

%bounds on input
param.opt.umax=50;
param.opt.umax1=50;
param.opt.umin1=-50;
param.opt.umin=-50;

param.opt.nX=4; %state size
param.opt.nU=1; %input size
% FILL IN cost matrices for the optimization problem
param.opt.Q  = diag([1, 1, 0.1, 0.1]);
param.opt.Qf = diag([10, 500, 1, 1]);
param.opt.R  = 0.01;

%time step bounds
param.opt.hmin=param.dt;
param.opt.hmax=param.dt;

end