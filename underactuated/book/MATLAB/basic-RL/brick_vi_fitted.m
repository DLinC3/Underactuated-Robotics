% ======================================================================
% Fitted Value Iteration with Neural Network Approximation
% Returns the optimal value function J and trained network
% ======================================================================
function J = brick_vi_fitted
global dt;

% Define state and action spaces
tic;
q_bins = -3:0.2:3;      % Discretization for position
qdot_bins = -4:0.2:4;   % Discretization for velocity
a = linspace(-1,1,9);   % Discrete action set
dt = 0.05;              % Time step
gamma = 0.99;           % Discount factor
converged = 2.5*1e-5;   % Convergence threshold

% Create the state grid
[q_mesh, qdot_mesh] = ndgrid(q_bins, qdot_bins);
s = [q_mesh(:)'; qdot_mesh(:)'];  % 2xN matrix
ns = size(s,2); 
na = numel(a);

% Generate all state-action pairs
S_all = repmat(s, 1, na);                   % Repeat states for each action
A_all = reshape(repelem(a, ns), 1, []);      % Repeat actions for each state

% Compute the next states (transition)
Sn = S_all + dynamics(S_all, A_all) .* dt;

% Initialize the neural network
layer_size1 = 20;
layer_size2 = 20;
layers = [ 
    featureInputLayer(2, 'Name', 'input')
    fullyConnectedLayer(layer_size1, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(layer_size2, 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(1, 'Name', 'output')
    regressionLayer('Name', 'regression')
];

% Training parameters
trainOpts = trainingOptions('adam', ...
    'MaxEpochs', 200, ...
    'InitialLearnRate', 0.01, ...
    'MiniBatchSize', 30, ...
    'Verbose', false);

C_all = cost(S_all, A_all);
C_reshaped = reshape(C_all, ns, na);
J_init = min(C_reshaped, [], 2);  % Initial target value

% Initial training to initialize the network
[trainedNet, info] = trainNetwork(s', J_init, layers, trainOpts);

% Main loop of value iteration
last_loss = inf;
figure(10);
for iter = 1:50
    % Three steps of the Bellman update
    J_next = predict(trainedNet, Sn');          % Step 1: Predict the next state's value
    J_next = J_next'; 

    C_all = cost(S_all, A_all);                  % Compute immediate cost
    Target = C_all + gamma .* J_next;            % Step 2: Compute the target value
    
    % Extract the minimum Q-value
    Target_matrix = reshape(Target, ns, na);
    Jnew = min(Target_matrix, [], 2);            % Step 3: Bellman optimality equation
    
    % Train the network
    [trainedNet, info] = trainNetwork(s', Jnew, layers, trainOpts);
    
    % Convergence check
    current_loss = info.TrainingLoss(end);
    err = abs(current_loss - last_loss);
    if err < converged
        break;
    end
    last_loss = current_loss;
    
    % Visualization
    J_pred = predict(trainedNet, s');
    vi_plot(J_pred, q_bins, qdot_bins, 10);
    disp(['Iter=', num2str(iter), ' Error=', num2str(err)]);
end

% Save the result
save('bricknet.mat', 'trainedNet');
J = predict(trainedNet, s');

% Simulation verification
simulate_control(trainedNet, q_bins, qdot_bins, gamma, dt);

end

% ==============================================================
% Control simulation function
% ==============================================================
function simulate_control(net, q_bins, qdot_bins, gamma, dt)
global dt;
T = 10; 
disp_dts = 5;
candidate_actions = linspace(-1,1,21);  % Fine action set

% Initialize trajectories
xtraj = zeros(2, T/dt);
utraj = zeros(1, T/dt);
x = clamp_state([-1; 1.2], q_bins, qdot_bins);  % Initial state

for i = 1:T/dt
    % Policy computation
    X = repmat(x, 1, numel(candidate_actions));
    U = candidate_actions;
    Sn = X + dynamics(X, U) .* dt;
    
    % Batch prediction of value function
    J_next = predict(net, Sn');
    C = cost(X, U);
    
    % Select the optimal action
    [~, idx] = min(C + gamma * J_next');
    u = candidate_actions(idx);
    
    % Update the state
    x = clamp_state(x + dynamics(x, u) .* dt, q_bins, qdot_bins);
    
    % Record trajectories
    xtraj(:, i) = x;
    utraj(i) = u;
    
    % Visualization
    if mod(i, disp_dts) == 0
        draw((i-1)*dt, x);
    end
end

% Plot the simulation results
plot_results(xtraj, utraj, T, dt);
end

% ==============================================================
% Helper functions
% ==============================================================
function xdot = dynamics(x, u)
    xdot = [x(2, :); u];
end

function C = cost(X, u)
    global dt;
    Q = diag([1, 1]) .* dt;
    R = 0.1 * dt;
    C = sum(X .* (Q * X), 1) + R * u.^2;
end

function s = clamp_state(s, q_bins, qdot_bins)
    q_lim = [min(q_bins), max(q_bins)];
    qdot_lim = [min(qdot_bins), max(qdot_bins)];
    
    s(1, :) = min(max(s(1, :), q_lim(1)), q_lim(2));
    s(2, :) = min(max(s(2, :), qdot_lim(1)), qdot_lim(2));
end

function vi_plot(J, q_bins, qdot_bins, fh)
    figure(fh);
    imagesc(q_bins, qdot_bins, reshape(J, numel(q_bins), [])');
    axis xy; colorbar;
    xlabel('Position'); ylabel('Velocity');
    title('Value Function');
    drawnow;
end

% ==============================================================
% Draw function for the brick (mass)
% ==============================================================
function draw(t, x)
persistent hFig blockx blocky;

if (isempty(hFig))
  hFig = figure(25);
  set(hFig, 'DoubleBuffer', 'on');
  blockx = [-1, -1, 1, 1, -1];
  blocky = [0, 0.5, 0.5, 0, 0];
end

figure(hFig);
clf;

% Draw the mass
brickcolor = [.75 .6 .5];
fill(blockx + repmat(x(1), 1, 5), blocky, brickcolor);
hold on;

faintline = [.6 .8 .65] * 1.1;
plot(min(blockx) + [0 0], [-5 5], 'k:', 'Color', faintline);
plot(max(blockx) + [0 0], [-5 5], 'k:', 'Color', faintline);

% Draw the ground
line([-5, 5], [0, 0], 'Color', [.3 .5 1], 'LineWidth', 1);
axis([-5 5 -1 2]);
axis equal;
title(['t = ', num2str(t)]);

drawnow;
end

function plot_results(xtraj, utraj, T, dt)
    t = 0:dt:T-dt;
    figure('Position', [100, 100, 800, 600]);
    
    subplot(3, 1, 1);
    plot(t, xtraj(1, :), 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Position');
    
    subplot(3, 1, 2);
    plot(t, xtraj(2, :), 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Velocity');
    
    subplot(3, 1, 3);
    plot(t, utraj, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Control');
    
    sgtitle('Control Trajectory');
end
