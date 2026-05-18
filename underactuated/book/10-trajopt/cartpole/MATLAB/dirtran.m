%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements direct transcription with SNOPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w, F, info] = dirtran(N, DYNAMICS, GRADIENTS, param)

nX = param.opt.nX;
nU = param.opt.nU;

%% Initialize decision variables
if isfield(param.opt, 'u0s')
    u0 = param.opt.u0s;
else
    u0 = zeros(nU*N, 1);
end
ulow  = param.opt.umin * ones(nU*N, 1);
uhigh = param.opt.umax * ones(nU*N, 1);

x0    = [];
xlow  = [];
xhigh = [];
for k = 1:nX
    x0    = [x0; linspace(param.opt.x0(k), param.opt.xd(k), N)];
    xlow  = [xlow; param.opt.xmin(k) * ones(1, N)];
    xhigh = [xhigh; param.opt.xmax(k) * ones(1, N)];
end
if isfield(param.opt, 'x0s')
    x0 = param.opt.x0s;
end

xlow(:, 1)    = param.opt.x0;
xhigh(:, 1)   = param.opt.x0;
xlow(:, end)  = param.opt.xfmin;
xhigh(:, end) = param.opt.xfmax;

hlow  = param.opt.hmin;
hhigh = param.opt.hmax;

w0    = [2/N; x0(:); u0];
wlow  = [hlow; xlow(:); ulow];
whigh = [hhigh; xhigh(:); uhigh];

%% Define constraint bounds
Flow  = [-inf; zeros((N-1)*nX, 1)];
Fhigh = [ inf; zeros((N-1)*nX, 1)];

idx_h = 1;
idxX  = @(k) (2 + (k-1)*nX) : (1 + k*nX);
idxU  = @(k) (2 + nX*N + (k-1)*nU) : (1 + nX*N + k*nU);

%% Build Jacobian sparsity pattern
iGfun = [];
jGvar = [];

num_vars = 1 + nX*N + nU*N;
iGfun = [iGfun; ones(num_vars, 1)];
jGvar = [jGvar; (1:num_vars)'];

for k = 1:(N-1)
    rowBlock = 1 + (k-1)*nX + (1:nX);

    iGfun = [iGfun; rowBlock'];
    jGvar = [jGvar; idxX(k+1)'];

    iGfun = [iGfun; repelem(rowBlock', nX)];
    jGvar = [jGvar; repmat(idxX(k)', nX, 1)];

    iGfun = [iGfun; repelem(rowBlock', nU)];
    jGvar = [jGvar; repmat(idxU(k)', nX, 1)];

    iGfun = [iGfun; rowBlock'];
    jGvar = [jGvar; repmat(idx_h, nX, 1)];
end

%% SNOPT parameters
setSNOPTParam('Major Iterations Limit', 10000);
setSNOPTParam('Minor Iterations Limit', 500);
setSNOPTParam('Major Optimality Tolerance', 1e-6);
setSNOPTParam('Major Feasibility Tolerance', 1e-6);
setSNOPTParam('Minor Feasibility Tolerance', 1e-6);
setSNOPTParam('Superbasics Limit', max(2000, num_vars));
setSNOPTParam('Derivative Option', 1);
setSNOPTParam('Verify Level', 0);
setSNOPTParam('Iterations Limit', 10000);

%% SNOPT inputs
wmul   = zeros(length(w0), 1);
wstate = zeros(length(w0), 1);
Fmul   = zeros(length(Flow), 1);
Fstate = zeros(length(Flow), 1);

USERFUN = @(w) userFun(w, N, DYNAMICS, GRADIENTS, param);

options.name      = 'Cartpole-DirTran';
options.start     = 'Cold';
options.screen    = 'on';
options.printfile = 'cartpole_dirtran_snopt.out';

A_struct = struct('row', [], 'col', [], 'val', []);
G_struct = struct('row', iGfun, 'col', jGvar);

ObjAdd = 0;
ObjRow = 1;

[w, F, info] = snopt(w0, wlow, whigh, wmul, wstate, ...
                     Flow, Fhigh, Fmul, Fstate, ...
                     USERFUN, ObjAdd, ObjRow, ...
                     A_struct, G_struct, options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for evaluating the cost, constraints, and gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, G] = userFun(w, N, DYNAMICS, GRADIENTS, param)

nX = param.opt.nX;
nU = param.opt.nU;

h = w(1);
X = reshape(w(2:(1+nX*N)), nX, N);
U = reshape(w(nX*N+2:end), nU, N);

f_dyn = zeros(nX, N-1);
J = 0;
Jsum_integrand = 0;
dJdx = zeros(nX, N);
dJdu = zeros(nU, N);

for i = 2:N
    xi = X(:, i-1);
    ui = U(:, i-1);

    xdot = DYNAMICS(0, xi, ui, param);
    f_dyn(:, i-1) = X(:, i) - (xi + h * xdot);

    Ji = cost(0, xi, ui, param);
    J = J + Ji * h;
    Jsum_integrand = Jsum_integrand + Ji;

    [dJdxi, dJdui] = cost_gradients(0, xi, ui, param);
    dJdx(:, i-1) = dJdxi.' * h;
    dJdu(:, i-1) = dJdui.' * h;
end

xN_error = X(:, N) - param.opt.xd;
J = J + final_cost(0, xN_error, [], param);

[dJdxN, ~] = final_cost_gradients(0, xN_error, [], param);
dJdx(:, N) = dJdxN.';
dJdu(:, N) = zeros(nU, 1);

dJdh = Jsum_integrand;

f = [J; f_dyn(:)];

values = [];
values = [values; dJdh];
values = [values; dJdx(:)];
values = [values; dJdu(:)];

for k = 1:(N-1)
    xi = X(:, k);
    ui = U(:, k);

    xdot = DYNAMICS(0, xi, ui, param);
    [dfdx, dfdu] = GRADIENTS(0, xi, ui, param);

    values = [values; ones(nX, 1)];

    block = -eye(nX) - h * dfdx;
    block_transposed = block.';
    values = [values; block_transposed(:)];

    dfdu_scaled = -h * dfdu;
    dfdu_transposed = dfdu_scaled.';
    values = [values; dfdu_transposed(:)];

    values = [values; -xdot];
end

G = values;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the final cost gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dJdx, dJdu] = final_cost_gradients(~, x, ~, param)
Qf = param.opt.Qf;
dJdx = 2 * x.' * Qf;
dJdu = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the final cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = final_cost(~, x, ~, param)
Qf = param.opt.Qf;
J = x.' * Qf * x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the running cost gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dJdx, dJdu] = cost_gradients(~, x, u, param)
Q = param.opt.Q;
R = param.opt.R;
dJdx = 2 * x.' * Q;
dJdu = 2 * u.' * R;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the running cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = cost(~, x, u, param)
Q = param.opt.Q;
R = param.opt.R;
J = x.' * Q * x + u.' * R * u;
end

function setSNOPTParam(paramstring, value)
snset([paramstring, '=', num2str(value)]);
end
