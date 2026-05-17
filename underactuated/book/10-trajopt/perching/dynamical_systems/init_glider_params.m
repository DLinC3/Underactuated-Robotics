%==========================================================================
%Contains all state estimation parameters
%==========================================================================
function param=init_param()

%x,z,theta,phi_dot,x_dot,z_doth,theta_dot
param.x0=[-3.4 0.1 0 0 7.4 0 0]';

param.dt=1/119; %time step
param.g=9.81;   %m/s
param.rho= 1.292; %m/s

dihedral = deg2rad(0); % dihedral angle
cm=0;
param.S_w = 0.05*cos(dihedral) + 0.035 + 0.0035; % wing + fus + tail
param.S_e = 0.0147;
param.l_w = 0.0+cm;
param.l_e = 0.022;
param.l_h = 0.27-cm;
param.m = 0.12;
param.I=0.0015;
 
param.cg_x=-0.065;
param.cg_z=-0.025;
param.vicon_offset=-0.052;

end
