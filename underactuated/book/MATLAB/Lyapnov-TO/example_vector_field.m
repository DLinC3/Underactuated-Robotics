function example_vector_field
%this MATLAB script demonstrates how to plot a vector field for a two dimensional dynamic system.
clear all;
close all;

%create vectors which define the range of your mesh points
x1=-4:0.1:4;
x2=-4:0.1:4;

%use the MATLAB command "meshgrid" to generate mesh points
[X1,X2] = meshgrid(x1,x2);

%set your input to zero
U=0;

%compute the derivatives of the dynamical system at each meshpoint using
%your dynamics function.
[DX1 DX2] = dynamics(X1,X2,U);

%normalize the derivatives
DX1=DX1./sqrt(DX1.^2+DX2.^2);
DX2=DX2./sqrt(DX1.^2+DX2.^2);

%use the MATLAB function "quiver" to plot the normalized vector field
figure(3)
hold on
quiver(X1,X2,DX1,DX2);
hold off

end

%system dynamics function (use pointwise multiplication)
function [Xdot1, Xdot2] = dynamics(X1,X2,U)


Xdot1=X2;
Xdot2=-X1+X1.*X2+U;

end


