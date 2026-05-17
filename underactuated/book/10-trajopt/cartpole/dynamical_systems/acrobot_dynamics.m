function xdot=acrobot_dynamics(t,x,u,param)

x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);

m1=param.m1;
m2=param.m2;
l1=param.l1;
lc1=param.lc1;
lc2=param.lc2;
b1=param.b1;
b2=param.b2;
I1=param.I1;
I2=param.I2;
w1=param.w1;
w2=param.w2;

g = 9.81;

xdot(1,:)=x(3);
xdot(2,:)=x(4);
% xdot(3,:) =(cos(x2)*sin(x2)*l1^2*lc2^2*m2^2*x3^2 + g*sin(x1 + x2)*cos(x2)*l1*lc2^2*m2^2 + I2*sin(x2)*l1*lc2*m2*x3^2 + 2*I2*sin(x2)*l1*lc2*m2*x3*x4 + I2*sin(x2)*l1*lc2*m2*x4^2 + b2*cos(x2)*l1*lc2*m2*x4 - u*cos(x2)*l1*lc2*m2 - I2*g*sin(x1)*l1*m2 - I2*b1*x3 + I2*b2*x4 - I2*u - I2*g*lc1*m1*sin(x1))/(- l1^2*lc2^2*m2^2*cos(x2)^2 + I2*l1^2*m2 + I1*I2)-w1*x3/sqrt(1+w2*x3^2);
% xdot(4,:) = -(I1*b2*x4 - I2*u - l1^2*m2*u - I2*b1*x3 - I1*u + I2*b2*x4 + b2*l1^2*m2*x4 + l1^3*lc2*m2^2*x3^2*sin(x2) + g*l1^2*lc2*m2^2*sin(x1 + x2) + I1*g*lc2*m2*sin(x1 + x2) - I2*g*l1*m2*sin(x1) - I2*g*lc1*m1*sin(x1) - 2*l1*lc2*m2*u*cos(x2) - b1*l1*lc2*m2*x3*cos(x2) + 2*b2*l1*lc2*m2*x4*cos(x2) + 2*l1^2*lc2^2*m2^2*x3^2*cos(x2)*sin(x2) + l1^2*lc2^2*m2^2*x4^2*cos(x2)*sin(x2) + g*l1*lc2^2*m2^2*sin(x1 + x2)*cos(x2) - g*l1^2*lc2*m2^2*cos(x2)*sin(x1) + I1*l1*lc2*m2*x3^2*sin(x2) + I2*l1*lc2*m2*x3^2*sin(x2) + I2*l1*lc2*m2*x4^2*sin(x2) + 2*l1^2*lc2^2*m2^2*x3*x4*cos(x2)*sin(x2) + 2*I2*l1*lc2*m2*x3*x4*sin(x2) - g*l1*lc1*lc2*m1*m2*cos(x2)*sin(x1))/(- l1^2*lc2^2*m2^2*cos(x2)^2 + I2*l1^2*m2 + I1*I2);

xdot(3,:) =(I2*b2*x4*(w2*x3^2 + 1)^(1/2) - I2*w1*x3 - I2*b1*x3*(w2*x3^2 + 1)^(1/2) - I2*u*(w2*x3^2 + 1)^(1/2) - I2*g*l1*m2*sin(x1)*(w2*x3^2 + 1)^(1/2) - I2*g*lc1*m1*sin(x1)*(w2*x3^2 + 1)^(1/2) - l1*lc2*m2*u*cos(x2)*(w2*x3^2 + 1)^(1/2) + I2*l1*lc2*m2*x3^2*sin(x2)*(w2*x3^2 + 1)^(1/2) + I2*l1*lc2*m2*x4^2*sin(x2)*(w2*x3^2 + 1)^(1/2) + b2*l1*lc2*m2*x4*cos(x2)*(w2*x3^2 + 1)^(1/2) + l1^2*lc2^2*m2^2*x3^2*cos(x2)*sin(x2)*(w2*x3^2 + 1)^(1/2) + g*l1*lc2^2*m2^2*sin(x1 + x2)*cos(x2)*(w2*x3^2 + 1)^(1/2) + 2*I2*l1*lc2*m2*x3*x4*sin(x2)*(w2*x3^2 + 1)^(1/2))/((w2*x3^2 + 1)^(1/2)*(- l1^2*lc2^2*m2^2*cos(x2)^2 + I2*l1^2*m2 + I1*I2));
 
 
xdot(4,:) =-(I1*b2*x4*(w2*x3^2 + 1)^(1/2) - I2*u*(w2*x3^2 + 1)^(1/2) - I2*w1*x3 - l1^2*m2*u*(w2*x3^2 + 1)^(1/2) - I2*b1*x3*(w2*x3^2 + 1)^(1/2) - I1*u*(w2*x3^2 + 1)^(1/2) + I2*b2*x4*(w2*x3^2 + 1)^(1/2) + b2*l1^2*m2*x4*(w2*x3^2 + 1)^(1/2) + I1*g*lc2*m2*sin(x1 + x2)*(w2*x3^2 + 1)^(1/2) - I2*g*l1*m2*sin(x1)*(w2*x3^2 + 1)^(1/2) - I2*g*lc1*m1*sin(x1)*(w2*x3^2 + 1)^(1/2) - 2*l1*lc2*m2*u*cos(x2)*(w2*x3^2 + 1)^(1/2) - l1*lc2*m2*w1*x3*cos(x2) + l1^3*lc2*m2^2*x3^2*sin(x2)*(w2*x3^2 + 1)^(1/2) + g*l1^2*lc2*m2^2*sin(x1 + x2)*(w2*x3^2 + 1)^(1/2) + I1*l1*lc2*m2*x3^2*sin(x2)*(w2*x3^2 + 1)^(1/2) + I2*l1*lc2*m2*x3^2*sin(x2)*(w2*x3^2 + 1)^(1/2) + I2*l1*lc2*m2*x4^2*sin(x2)*(w2*x3^2 + 1)^(1/2) - b1*l1*lc2*m2*x3*cos(x2)*(w2*x3^2 + 1)^(1/2) + 2*b2*l1*lc2*m2*x4*cos(x2)*(w2*x3^2 + 1)^(1/2) + 2*l1^2*lc2^2*m2^2*x3^2*cos(x2)*sin(x2)*(w2*x3^2 + 1)^(1/2) + l1^2*lc2^2*m2^2*x4^2*cos(x2)*sin(x2)*(w2*x3^2 + 1)^(1/2) + g*l1*lc2^2*m2^2*sin(x1 + x2)*cos(x2)*(w2*x3^2 + 1)^(1/2) - g*l1^2*lc2*m2^2*cos(x2)*sin(x1)*(w2*x3^2 + 1)^(1/2) + 2*I2*l1*lc2*m2*x3*x4*sin(x2)*(w2*x3^2 + 1)^(1/2) + 2*l1^2*lc2^2*m2^2*x3*x4*cos(x2)*sin(x2)*(w2*x3^2 + 1)^(1/2) - g*l1*lc1*lc2*m1*m2*cos(x2)*sin(x1)*(w2*x3^2 + 1)^(1/2))/((w2*x3^2 + 1)^(1/2)*(- l1^2*lc2^2*m2^2*cos(x2)^2 + I2*l1^2*m2 + I1*I2));
 
end