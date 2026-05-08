function [dfdx dfdu]=quadrotor_grads(t,x,u,w,param)

m=param.m;%mass
g=param.g;%gravity
L=param.L;%distance between quadrotor and cg
I=param.I;%inertia

dfdx=zeros(6,6);

% xdot(1,:) = x(4,:);
% xdot(2,:) = x(5,:);
% xdot(3,:) = x(6,:);
% xdot(4,:) = (-sin(x(3,:))/m)*(u(1)+u(2));
% xdot(5,:) = -g + (cos(x(3,:))/m)*(u(1)+u(2));
% xdot(6,:) = L/I*(-u(1)+u(2));


dfdx(1,4)=1;
dfdx(2,5)=1;
dfdx(3,6)=1;

dfdx(4,3)=(-cos(x(3))/m)*(u(1)+u(2));
dfdx(5,3)=(-sin(x(3))/m)*(u(1)+u(2));

dfdu=zeros(6,2);

dfdu(4,1)=(-sin(x(3))/m);
dfdu(5,1)=(cos(x(3))/m);
dfdu(6,1)=-L/I;

dfdu(4,2)=(-sin(x(3))/m);
dfdu(5,2)=(cos(x(3))/m);
dfdu(6,2)=L/I;

end