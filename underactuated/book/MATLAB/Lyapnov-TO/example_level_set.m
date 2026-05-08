function example_level_set
%this MATLAB script demonstrates how to generate a two dimensional plot of
%a level set of a function of the form V=x'Px, where x is a two dimensional
%vector x=[x1,x2] and V is the desired scalar level set value
clear all;
close all;

%define the center of the ellipse
x0=[0;0];

%define the first function's P matrix
P1=[12 -10;-10 12];

%define the second function's P matrix
P2=[10 -1;6 6];

%Define the desired level sets
V1=1;
V2=1;

%Plot the ellipses
figure(1)
hold on
plot_ellipse2D(P1,V1,x0,2,'r')
plot_ellipse2D(P2,V2,x0,2,'b')
axis([-1 1 -1 1])
hold off

end

function plot_ellipse2D(P,V,x0,linewidth,color)
%this MATLAB function plots a 2 Dimensional ellipse of the form
%x'Px-V=0.
%Inputs: P,V,x0,linewidth,color
% P: a 2x2 positive definite matrix
% V: a scalar level set value
% x0: a 2 dimensional vector specifying the origin of the ellipse
% linewidth: the desired line thickness of the ellipse
% color: the desired color of the ellipse
%Outputs:none

if min(eig(P))<0
    
    error('P is not positive definite')
    
end

%a vector of theta values
%we will use cylindrical components to plot our ellipse
THETA=0:.01:2*pi-0.01;

%the center of the ellipse
xc=x0(1);
zc=x0(2);

% the 2x2 P matrix
p11=P(1,1);
p22=P(2,2);
p12=P(1,2);
p21=P(2,1);


% for every theta from 0-2*pi, compute the radius of the ellipse
for k=1:length(THETA);
            
    theta=THETA(k);
    
    %compute the radius
    r= (V*p11*cos(theta)^2 + V*p22*sin(theta)^2 + V*p12*cos(theta)*sin(theta) + V*p21*cos(theta)*sin(theta))^(1/2)/(p11*cos(theta)^2 + p22*sin(theta)^2 + p12*cos(theta)*sin(theta) + p21*cos(theta)*sin(theta));
    
   %transform back to cartesian coordinates
   x(k)=r*cos(theta)+xc;
   z(k)=r*sin(theta)+zc;
       
    
end

%plot the ellipse
plot(x,z,color,'LineWidth',linewidth)

end

