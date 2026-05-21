% vanderpol_roa.m
%    ẋ₁ = -x₂
%    ẋ₂ = x₁ - (1 - x₁²)x₂
clear all; close all; clc;
%% (a)
figure(1);
plot_vector_field();
title('Vector Field');

%% (b)
% the linearization at 0 is
%   A = [0, -1; 1, -1].
A = [0 -1; 1 -1];
Q = eye(2);
P = lyap(A', Q)  % A'*P + P*A + Q = 0.
disp('eig(P) = ')
disp(eig(P))

%% (c)
% Create a grid
x_range = linspace(-4, 4, 3000); 
[X1, X2] = meshgrid(x_range, x_range);

%  V(x) = x' P x.
V = P(1,1)*X1.^2 + 2*P(1,2)*X1.*X2 + P(2,2)*X2.^2;

dV_dx1 = 2*P(1,1)*X1 + 2*P(1,2)*X2;
dV_dx2 = 2*P(1,2)*X1 + 2*P(2,2)*X2;

dx1 = -X2;
dx2 = X1 - (1 - X1.^2).*X2;

% Compute Vdot
dot_V = dV_dx1 .* dx1 + dV_dx2 .* dx2;

% Find points where dotV >= 0 and estimate the largest level set inside ROA
dot_V_non_neg = dot_V >= 0;
V_non_neg = V(dot_V_non_neg);

c_min_positive = min(V_non_neg);
disp('=====');
fprintf('Estimated maximum ρ from grid search: %.6f\n', c_min_positive);


% Plot contours: V = ρ (red) and dotV = 0 (blue).
figure(2);
plot_vector_field();  
hold on;
contour(X1, X2, V, [c_min_positive, c_min_positive], 'r', 'LineWidth', 2);
contour(X1, X2, dot_V, [0, 0], 'b', 'LineWidth', 2);
title('Grid-based ROA on Vector Field');
legend('Vector Field', 'V = ρ', 'dotV = 0', 'Location', 'best');
hold off;

%% (d)
% bilinear alternation
% The multiplier λ(x) is taken as x'*S*x

x=msspoly('x',2);
x1 = x(1);  x2 = x(2);

% using the same P
V_sym = P(1,1)*x1^2 + 2*P(1,2)*x1*x2 + P(2,2)*x2^2;

dVdx_sym = [diff(V_sym, x1); diff(V_sym, x2)];
f_sym = [-x2; x1 - (1 - x1^2)*x2];
Vdot_sym = dVdx_sym' * f_sym;

rho = msspoly('rho');
% Initialize ρ as fixed number 0.5
rho_fixed = 0.5;

max_iter = 30;
tol = 1e-5;

for iter = 1:max_iter
    % Step 1: For fixed ρ, find a multiplier λ(x) = x' S x

    s=msspoly('s',4);
    S = [ s(1) s(2) ; s(3) s(4) ];
    lambda_poly = x' * S * x;
    % slack variable gamma
    g = msspoly('g');

    constraint_slack = -Vdot_sym + lambda_poly*(V_sym - rho_fixed) - g;

    pr1 = mssprog;
    pr1.free = s;
    pr1.pos = g;

    pr1.sos = constraint_slack;
    pr1.sos = lambda_poly;
    pr1.sedumi = -g;
    % pr1({g});
    s_fixed = pr1({s});

    
    % Step 2: With λ fixed, maximize ρ 

    % S_fixed
    S_fixed = [ s_fixed(1) s_fixed(2) ; s_fixed(3) s_fixed(4) ];
    lambda_fixed = x' * S_fixed * x;
    pr2 = mssprog;
    pr2.pos = rho;
    pr2.sos = -Vdot_sym + lambda_fixed*(V_sym - rho);
    pr2.sedumi = -rho;
    rho_new = pr2({rho});
    
    fprintf('Iteration %d: rho = %.6f\n', iter, rho_new);
    if abs(rho_new - rho_fixed) < tol
        break;
    else
        rho_fixed = rho_new;
    end
end

disp('=====');
fprintf('SOS estimated maximum rho: %.6f\n', rho_new);

% Plot the SOS-based ROA ellipse on the vector field.
figure(3);
plot_vector_field();
hold on;
plot_ellipse2D(P, rho_new, [0;0], 2, 'g');
title('SOS-based ROA on Vector Field');
legend('Vector Field', 'SOS ROA', 'Location', 'best');
hold off;


figure(4);
plot_vector_field(); 
hold on;
plot_ellipse2D(P, rho_new, [0;0], 2, 'g');
title('Inside and Outside Trajecytory on Vector Field');

for j=1:10 % 10 inside trajectory
    x=linspace(-.71,.71,2001);
    for i=1:length(x);
        y1(i)=max(roots([1 -x(i) 1.5*x(i)^2-rho_new]));
        y2(i)=min(roots([1 -x(i) 1.5*x(i)^2-rho_new]));
    end
    
    idx=ceil(rand*2001);
    x0=[x(idx) y1(idx)];
    [t,z]=ode45(@(t,z) [-z(2); z(1)+z(2)*(z(1)^2-1)], [0 10], x0);
    plot(z(:,1),z(:,2),'k','LineWidth',.7);
end

c = 7;  % Outside trajectory as V(x) = 7
P = [1.5 -0.5; -0.5 1];

x = linspace(-3, 3, 2001);  
valid_points = [];
for i = 1:length(x)
    coeff = [1, -x(i), 1.5*x(i)^2 - c];
    roots_result = roots(coeff);
    real_roots = real(roots_result(abs(imag(roots_result)) < 1e-6));
    if ~isempty(real_roots)
        valid_points = [valid_points; x(i), max(real_roots); x(i), min(real_roots)];
    end
end

hold on;

options = odeset(...
    'RelTol', 1e-4,...
    'AbsTol', 1e-6,...
    'Events', @escapeEvent,...
    'MaxStep', 0.1);  % maxstep

for j = 1:10
    if isempty(valid_points)
        error('No valid initial points!');
    end
    idx = randi(size(valid_points,1));
    x0 = valid_points(idx,:);
    
    try
        [t,z] = ode15s(@vdp_reversed, [0 20], x0, options);
        
        if ~isempty(t)
            plot(z(:,1), z(:,2), 'Color', [0.9 0.2 0.2], 'LineWidth', 1.5);
            plot(x0(1), x0(2), 'o', 'MarkerFaceColor', [1 0.8 0], 'MarkerSize', 6);
        end
    catch ME
        fprintf('Trajectory %d failed: %s\n', j, ME.message);
    end
end

hold off;
axis([-4 4 -4 4]);  % fixed frame range
grid on;


%% Help Functions


function plot_vector_field()
    clf;
    hold on;
    
    x1 = -4:0.1:4;
    x2 = -4:0.1:4;
    [X1, X2] = meshgrid(x1, x2);
    
    [DX1, DX2] = dynamics(X1, X2, 0);
    
    mag = sqrt(DX1.^2 + DX2.^2);
    DX1 = DX1 ./ mag;
    DX2 = DX2 ./ mag;
    
    quiver(X1, X2, DX1, DX2);
    axis equal;
    xlabel('x₁'); 
    ylabel('x₂');
    grid on;
    hold off;
end

function plot_ellipse2D(P,V,x0,linewidth,color)
    if min(eig(P)) < 0
        error('P is not positive definite');
    end
    
    p11 = P(1,1);
    p22 = P(2,2);
    p12 = P(1,2);
    p21 = P(2,1);
    
    THETA = 0:.01:2*pi-0.01;
    
    x = zeros(size(THETA));
    z = zeros(size(THETA));
    
    xc = x0(1);
    zc = x0(2);
    
    for k = 1:length(THETA)
        theta = THETA(k);
        
        numerator = V*(p11*cos(theta)^2 + p22*sin(theta)^2 + ...
                      p12*cos(theta)*sin(theta) + p21*cos(theta)*sin(theta));
        denominator = p11*cos(theta)^2 + p22*sin(theta)^2 + ...
                     p12*cos(theta)*sin(theta) + p21*cos(theta)*sin(theta);
        
        r = sqrt(numerator)/denominator;
        
        x(k) = r*cos(theta) + xc;
        z(k) = r*sin(theta) + zc;
    end
    
    plot(x, z, color, 'LineWidth', linewidth);
end

function [Xdot1, Xdot2] = dynamics(X1, X2, ~)
    Xdot1 = -X2;
    Xdot2 = X1 - (1 - X1.^2).*X2;
end

% Inverse program definition (separate function) (for outside trajectory)
function dzdt = vdp_reversed(t,z)
    dzdt = [-z(2); 
            z(1) + z(2)*(z(1)^2 - 1)];
end

function [value,isterminal,direction] = escapeEvent(t,z)
    % since the plot is between x1,x2=[-4,4]
    value = max(abs(z)) - 5;  % Stop when trajectory exceeds 5 units
    isterminal = 1;
    direction = 1;
end
