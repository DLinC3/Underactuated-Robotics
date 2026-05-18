% originally from AEM 5431 notes
clear; clc; close all;

param = struct; % structure variable for all parameters
param.nx = 4; % # of states
param.nu = 2; % # of inputs
param.x0hat = zeros(param.nx, 1); % initial state
param.xfhat = [10; 10; 0; 0]; % final state
param.tmin = 5; % min time
param.tmax = 20; % max time
param.N = 10; % # of discretization points
param.taustamps = linspace(0, 1, param.N); % time stamps
param.umax = 1*ones(param.nu, 1); % upper bound on input
param.umin = -param.umax; % lower bound on input
param.velmax = 3; % maximum speed
param.obst_cen = [5; 4.5]; % center location of the obstacle
param.obst_rad = 3; % radius of the obstacle
% parameters of the continuous-time dynamics
param.Ac = [zeros(2, 2), eye(2); zeros(2, 4)];
param.Bc = [zeros(2, 2); eye(2)];
param.xbar2pos = [eye(2), zeros(2, 2), zeros(2, 2)]; % matrix that maps the augmented state to position
param.xbar2vel = [zeros(2, 2), eye(2), zeros(2, 2)]; % matrix that maps the augmented state to velocity
Ak = zeros(param.nx+2, param.nx+2, param.N-1); % parameters for linearized dynamics: Ak
Bk = zeros(param.nx+2, param.nu+1, param.N-1); % parameters for linearized dynamics: Bk
ck = zeros(param.nx+2, param.N-1); % parameters for linearization dynamics: ck
%%
maxiter = 300; % max iters for prox-linear
eps = 5e-3; % stopping tolerance for prox-linear
gam = 1e2; % exact penalty parameter
rho = 1/gam; % weight parameter in prox-linear
eps_path = 1e-6; % relaxation parameters for path constraint violation
xi_tilde = rand(param.nx+2, param.N);
% linearization point of the state sequence
eta_tilde = rand(param.nu+1, param.N-1);
% linearization point of the input sequence
for iter = 1:maxiter
    % linearization
    for k = 1:param.N-1
        init_val = [eta_tilde(:, k); xi_tilde(:, k); reshape(eye(param.nx+2), [], 1); zeros((param.nx+2)*(param.nu+1), 1)];
        [~, int_hist] = ode45(@(t, y) mt_dyn_ode(t, y, param), [param.taustamps(k), param.taustamps(k+1)], init_val);
        coeff = int_hist(end, :)';
        xbar = coeff(param.nu+2:param.nu+param.nx+3);
        Ak(:, :, k) = reshape(coeff(param.nu+param.nx+4:param.nu+1+(param.nx+2)*(param.nx+3)), param.nx+2, param.nx+2);
        Bk(:, :, k) = reshape(coeff(param.nu+2+(param.nx+2)*(param.nx+3):param.nu+1+(param.nx+2)*(param.nu+param.nx+4)), param.nx+2, param.nu+1);
        ck(:, k) = -Ak(:, :, k)*xi_tilde(:, k) - Bk(:, :, k)*eta_tilde(:, k) + xbar;
    end
    % optimizing via prox-linear
    yalmip('clear')
    % variables for trajectory optimization
    xi = sdpvar(param.nx+2, param.N, 'full');
    % variables for the state trajectory
    eta = sdpvar(param.nu+1, param.N-1, 'full');
    % variables for the input trajectory
    q = sdpvar(param.nx+2, param.N-1, 'full');
    % slack variables for relaxing dynamics constraints
    r = sdpvar(param.nx+2, param.N-1, 'full');
    % slack variables for relaxing dynamics constraints
    % dynamics constraints for trajectory optimization
    constr = [xi(:, 1) == [param.x0hat; zeros(2, 1)], ...
              xi(1:param.nx, param.N) == param.xfhat];
    % boundary condition
    for k = 1:param.N-1
        constr = [constr, xi(:, k+1) == Ak(:, :, k)*xi(:, k) + Bk(:, :, k)*eta(:, k) + ck(:, k) + q(:, k) - r(:, k)];
        % dynamics constraints
        constr = [constr, eta(:, k) <= [param.umax; param.tmax], ...
                  eta(:, k) >= [param.umin; param.tmin], ...
                  xi(param.nx+1, k+1) <= eps_path];
        % path constraints
        constr = [constr, norm(param.xbar2vel*xi(:, k+1)) <= param.velmax];
    end
    constr = [constr, q >= zeros(param.nx+2, param.N-1), r >= zeros(param.nx+2, param.N-1)];
    obj = xi(end, end) + gam*sum(sum(q+r)) + (0.5/rho)*(sum(sum((xi-xi_tilde).^2)) + sum(sum((eta-eta_tilde).^2)));
    options = sdpsettings('verbose', 0, 'solver', 'mosek');
    solution = optimize(constr, obj, options);
    xi_opt = value(xi); % optimal state trajectory
    eta_opt = value(eta); % optimal input trajectory
    constr_vio = sum(sum(value(q) + value(r))); % compute constraint violation
    step = (sum(sum((xi_opt-xi_tilde).^2)) + sum(sum((eta_opt-eta_tilde).^2)))^0.5; % compute step length
    if step <= eps && constr_vio <= eps
        break
    else
        xi_tilde = xi_opt;
        eta_tilde = eta_opt;
        fprintf('iteration %d, slack magnitude %f, step length %f \n', iter, constr_vio, step);
    end
end

%%
xopt = zeros(param.nx, param.N);
xopt(:, 1) = param.x0hat;
x_hist = [];
t_hist = [];
ts = zeros(1, param.N);
for k = 1:param.N-1
    ts(k+1) = ts(k) + eta_opt(end, k)*(param.taustamps(k+1)-param.taustamps(k));
end

for k = 1:param.N-1
    [t, int_hist] = ode45(@(t, x) reshape(param.Ac*x + param.Bc*eta_opt(1:param.nu, k), [], 1), [ts(k), ts(k+1)], xopt(:, k));
    x_hist = [x_hist, int_hist(2:end, :)'];
    t_hist = [t_hist, t(2:end, :)'];
    xopt(:, k+1) = x_hist(:, end);
end

figure('Position', [100, 100, 700, 700])
circles(param.obst_cen(1), param.obst_cen(2), param.obst_rad, 'edgecolor', 0.8*ones(1, 3), 'facecolor', 0.8*ones(1, 3), 'LineWidth', 1)
hold on
plot(x_hist(1, :), x_hist(2, :), 'b', 'LineWidth', 2)
plot(xopt(1, :), xopt(2, :), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', [1 .6 .6])
hold off
vel_hist = x_hist(3:4, :);
figure('Position', [100, 100, 700, 700])
plot(t_hist, sum(vel_hist.^2).^0.5, 'b', 'LineWidth', 2)
hold on
plot(ts, sum(xopt(3:4, :).^2).^0.5, 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', [1 .6 .6])
ylim([0 1.5*param.velmax])

function dydt = mt_dyn_ode(t, y, param)
u = y(1:param.nu); % extract input
s = y(param.nu+1);
xbar = y(param.nu+2:param.nu+param.nx+3);
% extract augmented state
x = xbar(1:param.nx); % extract state
Phi_x = reshape(y(param.nu+param.nx+4:param.nu+1+(param.nx+2)*(param.nx+3)), param.nx+2, param.nx+2);
% extract state-to-state Jacobian
Phi_u = reshape(y(param.nu+2+(param.nx+2)*(param.nx+3):param.nu+1+(param.nx+2)*(param.nu+param.nx+4)), param.nx+2, param.nu+1);
% extract input-to-state Jacobian
% parameters for ODEs
PfbarPxbar = s*[param.Ac, zeros(param.nx, 2); ...
    ((1/param.velmax^2)*norm(param.xbar2vel*xbar)^2 - 1 > 0)*(1/param.velmax^2)*(param.xbar2vel*xbar)'*param.xbar2vel + ...
    (1 - (1/param.obst_rad^2)*norm(param.xbar2pos*xbar - param.obst_cen)^2 > 0)*(-(1/param.obst_rad^2)*(param.xbar2pos*xbar - param.obst_cen)'*param.xbar2pos); ...
    zeros(1, param.nx+2)];
% continuous-time derivatives: partial fbar / partial xbar
PfbarPubar = [s*param.Bc, param.Ac*x + param.Bc*u; ...
    zeros(1, param.nu), 0.5*max((1/param.velmax^2)*norm(param.xbar2vel*xbar)^2 - 1, 0)^2 + ...
    0.5*max(1 - (1/param.obst_rad^2)*norm(param.xbar2pos*xbar - param.obst_cen)^2, 0)^2; ...
    zeros(1, param.nu), 1];
% continuous-time derivatives: partial fbar / partial u
% ODEs
DuDt = zeros(param.nu+1, 1); % ZOH: constant input
DxbarDt = s*[param.Ac*x + param.Bc*u; ...
    0.5*max((1/param.velmax^2)*norm(param.xbar2vel*xbar)^2 - 1, 0)^2 + ...
    0.5*max(1 - (1/param.obst_rad^2)*norm(param.xbar2pos*xbar - param.obst_cen)^2, 0)^2; ...
    1]; % continuous-time dynamics for augmented state
DPhi_xDt = PfbarPxbar*Phi_x;
% continuous-time dynamics for Jacobians: Ak
DPhi_uDt = PfbarPxbar*Phi_u + PfbarPubar;
% continuous-time dynamics for Jacobians: Bk
% reshape into vector
dydt = [DuDt; DxbarDt; reshape(DPhi_xDt, [], 1); reshape(DPhi_uDt, [], 1)]; % stack all derivatives into one vector
end

function h = circles(x, y, r, varargin)
%CIRCLES Draw one or more filled circles using PATCH.
%   circles(x, y, r, ...) draws circles centered at (x, y) with radius r.
%   Extra name-value arguments are forwarded to PATCH, e.g.
%   circles(0, 0, 1, 'edgecolor', [0.8 0.8 0.8], 'facecolor', [0.8 0.8 0.8]).

theta = linspace(0, 2*pi, 200);
x = x(:);
y = y(:);
r = r(:);
n = max([numel(x), numel(y), numel(r)]);

if isscalar(x)
    x = repmat(x, n, 1);
end
if isscalar(y)
    y = repmat(y, n, 1);
end
if isscalar(r)
    r = repmat(r, n, 1);
end

h = gobjects(n, 1);
for i = 1:n
    xx = x(i) + r(i)*cos(theta);
    yy = y(i) + r(i)*sin(theta);
    h(i) = patch(xx, yy, 'w', varargin{:});
end
axis equal

if nargout == 0
    clear h
end
end
