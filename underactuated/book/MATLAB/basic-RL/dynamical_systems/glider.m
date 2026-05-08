function xdot = glider(t,x,u,param)
%global m g rho S_w S_e I l_h l_w l_e;

%u=ppval(t,utape);

Sw=param.S_w;
Se=param.S_e;
l=param.l_h;
lw=param.l_w;
le=param.l_e;

m=param.m;
g=param.g;
rho=param.rho;
I=param.I;

u(2)=0;

 q1 = x(1,:); q2 = x(2,:);  q3 = x(3,:); q4 = x(4,:);
 qdot1 = x(5,:); qdot2 = x(6,:); qdot3 = x(7,:); qdot4=u(1);

xwdot = qdot1 - lw*qdot3.*sin(q3);
zwdot=  qdot2 + lw*qdot3.*cos(q3);
alpha_w = q3 - atan2(zwdot,xwdot);
Fw = rho*Sw*sin(alpha_w).*(zwdot.^2+xwdot.^2);

xedot = qdot1 + l*qdot3.*sin(q3) + le*(qdot3+qdot4).*sin(q3+q4); 
zedot = qdot2 - l*qdot3.*cos(q3) - le*(qdot3+qdot4).*cos(q3+q4);
alpha_e = q3+q4-atan2(zedot,xedot);
Fe = rho*Se*sin(alpha_e).*(zedot.^2+xedot.^2);

xdot=x;

xdot(4,:)=u(1);

xdot(5,:)=(-(Fw.*sin(q3) + Fe.*sin(q3+q4)+u(2)*cos(q3))./m);

xdot(6,:)=((Fw.*cos(q3) + Fe.*cos(q3+q4)+u(2)*sin(q3))./m -g);

xdot(7,:)=((Fw.*lw - Fe.*(l*cos(q4)+le))./I);

xdot(6,:)=((Fw.*cos(q3) + Fe.*cos(q3+q4))./m-g)+u(2)*sin(q3)./m;

xdot(1,:)=x(5,:);

xdot(2,:)=x(6,:);

xdot(3,:)=x(7,:);

end