function J = pend_vi_fitted
global dt;

tic;
% Define the mesh for theta and theta_dot.
q_bins = linspace(0, 2*pi, 25);   
qdot_bins = linspace(-10, 10, 25);
if (isempty(q_bins) || isempty(qdot_bins))
    error('you need to define the mesh');
end

% Define the discrete actions.
a = linspace(-2, 2, 10);
if (isempty(a))
    error('you need to define the action set');
end

% Set dynamics time-step.
dt = 0.05;

% Create the mesh.
[q, qdot] = ndgrid(q_bins, qdot_bins);
s = [reshape(q, 1, []); reshape(qdot, 1, [])]; % 2 x ns
ns = size(s, 2);
na = numel(a);

% Generate all state-action pairs.
S = repmat(s, 1, na);               % each state repeated for every action
A = reshape(repmat(a, ns, 1), 1, []); % each action for every state

% Compute one-step dynamics.
Sn = S + dynamics(S, A).*dt;
% Apply angle wrapping on the first state (angle).
Sn(1,:) = mod(Sn(1,:), 2*pi);

% Setup value iteration.
J = zeros(ns, 1);  % arbitrary initialization
gamma = 0.99;      % discount factor
loadFromData = 0;  % set to 1 if you want to load a pretrained network

fh = figure(10);

if(loadFromData == 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the neural network architecture
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    layer_size1 = 64;
    layer_size2 = 64;
    layers = [ 
        featureInputLayer(2, 'Name','input')
        fullyConnectedLayer(layer_size1, 'Name','fc1')
        reluLayer('Name','relu1')
        fullyConnectedLayer(layer_size2, 'Name','fc2')
        reluLayer('Name','relu2')
        fullyConnectedLayer(1, 'Name','output')
        regressionLayer('Name','regression')
    ];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set training options for the network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trainOpts = trainingOptions('adam', ...
        'MaxEpochs', 300, ...
        'InitialLearnRate', 0.01, ...
        'MiniBatchSize', 30, ...
        'Verbose', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the network with initial target values.
    % Compute one–step cost for every (state,action) pair.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_all = cost(S, A);
    C_matrix = reshape(C_all, ns, na);
    J_init = min(C_matrix, [], 2);  % take the minimum cost over actions
    [trainedNet, info] = trainNetwork(s', J_init, layers, trainOpts);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Value iteration loop with neural network update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    last_loss = inf;
    iter = 1;
    Jlast = J;
    max_iter = 200;
    while (iter <= max_iter)
        % 1) Evaluate the value function at the next state (x_{k+1})
        J_next = predict(trainedNet, Sn');
        J_next = reshape(J_next, ns, na);
        
        % 2) Apply the Bellman update: compute target Q-values.
        Target_matrix = C_matrix + gamma .* J_next;
        Jnew = min(Target_matrix, [], 2);  % select the best (lowest) cost

        % 3) Update the approximate value function by retraining.
        [trainedNet, info] = trainNetwork(s', Jnew, layers, trainOpts);
        
        % Compute convergence metrics.
        loss = info.TrainingLoss(end);
        err = max(abs(Jnew - Jlast));
        disp(['iteration = ', num2str(iter), ' ; training loss diff = ', num2str(abs(loss - last_loss)), ' ; max change in J = ', num2str(err)]);
        if err < 1e-5
            break;
        end
        last_loss = loss;
        Jlast = Jnew;
        J = Jnew;
        vi_plot(J, q_bins, qdot_bins, fh);
        iter = iter + 1;
    end

    disp('Value Estimate converged!');
    xlabel('q'); ylabel('q_{dot}'); title('Value function'); colorbar;
    set(gcf, 'PaperPosition', [0 0 8 5]);
    set(gcf, 'PaperSize', [8 5]);
    saveas(gcf, 'mintimevi_pend', 'pdf');
    save('pendnet.mat', 'trainedNet');
else
    temp = load('pendnet.mat', 'trainedNet');
    trainedNet = temp.trainedNet;
end

% Simulation
toc;
disp('Press Enter to simulate...'); pause;
T = 10;
disp_dts = 1;
candidate_actions = linspace(-2,2,21);
for j = 1:1
    % Start from a random initial state (normalized).
    x = normalize(randn(2,1), q_bins, qdot_bins);
    xtraj = zeros(2, T/dt);
    for i = 1:(T/dt)
        xtraj(:,i) = x;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the policy using the approximate value function.
        % For the current state x, evaluate candidate actions.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_candidates = numel(candidate_actions);
        X = repmat(x, 1, num_candidates);
        Sn_candidate = X + dynamics(X, candidate_actions).*dt;
        Sn_candidate(1,:) = mod(Sn_candidate(1,:), 2*pi);
        J_candidate = predict(trainedNet, Sn_candidate');
        J_candidate = J_candidate(:);
        C_candidate = cost(X, candidate_actions);
        Q_candidate = C_candidate + gamma * J_candidate;
        [~, idx] = min(Q_candidate);
        u = candidate_actions(idx);
        
        if mod(i, disp_dts) == 0
            draw((i-1)*dt, x);
        end
        
        % Update state using dynamics.
        x = normalize(x + dynamics(x, u).*dt, q_bins, qdot_bins);
    end
end

end  % end of pend_vi_fitted

% =============================================================
% Continuous dynamics for the pendulum.
% m = 1, l = 0.5, b = 0.1, I = 0.25, g = 9.8.
% =============================================================
function xdot = dynamics(x, u)
    m = 1;   
    l = 0.5;  
    b = 0.1; 
    I = 0.25;
    g = 9.8;
    xdot = [x(2,:); (u - m*g*l*sin(x(1,:)) - b*x(2,:)) / I];
end

% ===================================================================
% Instantaneous cost function.
% We penalize the angular error from upright (pi), angular velocity, and control effort.
% ===================================================================
function C = cost(X, u)
global dt
[Q, R] = get_QR;
% Compute angular error so that error is in [-pi,pi].
theta_error = mod(X(1,:) - pi + pi, 2*pi) - pi;
C = (Q(1,1) * (theta_error).^2 + Q(2,2) * (X(2,:)).^2) + R * (u).^2;
end

function [Q, R] = get_QR
global dt;
Q = diag([10 10]) * dt;
R = 1 * dt;
end

% ==============================================================
% Normalizes the state by wrapping the angle and clamping velocity.
% ==============================================================
function s = normalize(s, q_bins, qdot_bins)
N = size(s,2);
s(1,:) = mod(s(1,:), 2*pi);
smax = repmat([q_bins(end); qdot_bins(end)], 1, N);
smin = repmat([q_bins(1); qdot_bins(1)], 1, N);
ind = s > smax;
s(ind) = smax(ind);
ind = s < smin;
s(ind) = smin(ind);
end

% ===============================================================
% Plot the value function.
% q_bins and qdot_bins are vectors and fh is the figure handle.
% ===============================================================
function vi_plot(J, q_bins, qdot_bins, fh)
set(0, 'CurrentFigure', fh);
clf;
n1 = numel(q_bins);
n2 = numel(qdot_bins);
imagesc(q_bins, qdot_bins, reshape(J, n1, n2)');
axis xy; colorbar;
xlabel('q'); ylabel('q_{dot}');
title('Value Function');
drawnow;
end

% ==============================================================
% Draw function for the pendulum.
% ==============================================================
function draw(t, x)
persistent hFig base a1 a2 ac1 ac2 raarm;
if isempty(hFig)
    hFig = figure(25);
    set(hFig, 'DoubleBuffer', 'on');
    
    % parameters for drawing (feel free to adjust for aesthetics)
    a1 = 0.75;  
    ac1 = 0.415;
    av = pi*[0:.05:1];
    rb = 0.03; hb = 0.07;
    aw = 0.01;
    base = rb*[1 cos(av) -1 1; -hb/rb sin(av) -hb/rb -hb/rb]';
    arm = [aw*cos(av-pi/2) -a1+aw*cos(av+pi/2)
           aw*sin(av-pi/2) aw*sin(av+pi/2)]';
    raarm = [(arm(:,1).^2 + arm(:,2).^2).^.5, atan2(arm(:,2), arm(:,1))];
end

figure(hFig); cla; hold on; view(0,90);
patch(base(:,1), base(:,2), ones(size(base,1),1), 'b', 'FaceColor', [.3 .6 .4]);
patch(raarm(:,1).*sin(raarm(:,2)+x(1)-pi), ...
      -raarm(:,1).*cos(raarm(:,2)+x(1)-pi), zeros(size(raarm,1),1), 'r', 'FaceColor', [.9 .1 0]);
plot3(ac1*sin(x(1)), -ac1*cos(x(1)), 1, 'ko', 'MarkerSize',10,'MarkerFaceColor','b');
plot3(0, 0, 1.5, 'k.');
title(['t = ', num2str(t, '%.2f'), ' sec']);
set(gca, 'XTick',[], 'YTick',[]);
axis image; axis([-1.0 1.0 -1.0 1.0]);
drawnow;
end
