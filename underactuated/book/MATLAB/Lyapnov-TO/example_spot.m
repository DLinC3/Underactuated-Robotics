function example_spot
%This MATLAB script demonstrates how to use SPOT and SeDuMi to find the
%smallest value of "r" such that the equation
%4*(x^4)*(y^6)+r*(x^2)-x*(y^2)+y^2 is still a sum of squares.

%define your mss variables
x=msspoly('x');
y=msspoly('y');
r=msspoly('r');

%formulate the equation
q=4*(x^4)*(y^6)+r*(x^2)-x*(y^2)+y^2;

%initialize a mss program
pr=mssprog;

%set r, the variable you are searching for as a free variable
pr.free=r;

%constrain q to be a sum of squares
pr.sos=q;

%tell spot to minimize r
pr.sedumi=r;

%display r
pr({r})



end

