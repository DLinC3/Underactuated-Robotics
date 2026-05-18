function xdot=quadrotor(t,x,u,w,param)

m=param.m;%mass
g=param.g;%gravity
L=param.L;%distance between quadrotor and cg
I=param.I;%inertia

%Compute forward dynamics
xdot(1,:) = x(4,:);
xdot(2,:) = x(5,:);
xdot(3,:) = x(6,:);
xdot(4,:) = (-sin(x(3,:))/m)*(u(1)+u(2));
xdot(5,:) = -g + (cos(x(3,:))/m)*(u(1)+u(2));
xdot(6,:) = L/I*(-u(1)+u(2));

end

