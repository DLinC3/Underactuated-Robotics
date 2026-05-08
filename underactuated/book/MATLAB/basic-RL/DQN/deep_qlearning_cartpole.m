function deep_qlearning_cartpole
    % Global time step
    global dt;
    dt = 0.1;
    
    % [x; theta; xdot; thetadot]
    xd = [0; pi; 0; 0];
    
    % Discount factor
    gamma = 0.8;
    
    Qm.a = linspace(-2,2,10);
    
    %% Define the neural network for the Q function.
    % input dimension from 2 to 4.
    layers = [ ...
        featureInputLayer(4, "Name", "myFeatureInputLayer")
        fullyConnectedLayer(128, "Name", "myFullyConnectedLayer1")
        reluLayer("Name", "myReLu1")
        fullyConnectedLayer(64, "Name", "myFullyConnectedLayer2")
        reluLayer("Name", "myReLu2")
        fullyConnectedLayer(length(Qm.a), "Name", "myFullyConnectedLayer3") % 10 outputs
        regressionLayer("Name", "myRegressionLayer")
    ];
    %% Define options for the training
    opts = trainingOptions('adam', ...
        'MaxEpochs',100, ...
        'InitialLearnRate',0.01, ...
        'MiniBatchSize',1024, ...
        'Verbose',false, ...
        Plots="none");
    % create two networks, one for training and one for evaluation (target)
    Qm.layers = layers;% main network for training (main network)
    Qm.opts = opts;
    

    Qt = Qm;%network for evaluation (target network)
    
    x_fixed = 0; xdot_fixed = 0;
    theta_bins = linspace(0,2*pi,25);
    thetadot_bins = linspace(-10,10,25);
    [theta, thetadot] = ndgrid(theta_bins, thetadot_bins);
    s = [repmat(x_fixed, 1, numel(theta)); ...
         reshape(theta, 1, []); ...
         repmat(xdot_fixed, 1, numel(theta)); ...
         reshape(thetadot, 1, [])];
     
    x_bins = linspace(-2.5,2.5,25);
    xdot_bins = linspace(-10,10,25);
    thetadot_bins = linspace(-10,10,25);
    
    N = 20;      %train network every N simulation steps
    Ninit = 1000;  %number of simulation steps before training starts
    Nepisode = 200; % maximum steps per episode
    
    x0k = zeros(4, N);
    a0k = zeros(1, N);
    Qmk = zeros(N, length(Qm.a));
    Qm = learn_q_factor(x0k, Qmk, Qm);
    Qt = learn_q_factor(x0k, Qmk, Qm);
    
    fh = figure(10);
    Jlast = zeros(size(s,2), 1);
    dJnorm = [];
    kJ = [];
    
    % from ilqr code
    param.mc = 10;
    param.mp = 2;
    param.l = 0.5;
    param.g = 9.8;
    param.b = 0.1;
    param.d = 0.1;
    
    loadFromData = false;
    if ~loadFromData
        alpha = 0.8;%learning rate
        M = 400;%number of episodes
        P = 100;%number of steps before copying weights to target network
        k = 1;
        for i = 1:M

            x0 = [0; pi; 0; 0];


            done = false;
 
           epsilon = 0.01 + (1 - 0.01) * exp(-0.01 * i);
            j = 1;
            % while the episode is not complete
            
            while ~done
                display('epsiode || steps in episode || simulation steps')
                [i j k]
                
                draw_cartpole(0, x0, param);

                x0k(:,k) = x0;
                
                a_ind = get_action(Qt, x0, epsilon);
                a0k(:,k) = a_ind;
                
                x1 = x0 + cartpole_dynamics(x0, Qt.a(a_ind), param) * dt;
                
                if (abs(x1(1)) > max(abs(x_bins))) || (abs(x1(3)) > max(abs(xdot_bins))) || (j > Nepisode)
                    done = true;
                end
                j = j + 1;
                
                x1 = normalize(x1, x_bins, theta_bins, xdot_bins, thetadot_bins);
                x1k(:,k) = x1;
                % set
                x0 = x1;
                
                %train the main network every N simulation steps
                if mod(k, N) == 0 && k > Ninit
                    Qm = train(x0k, a0k, x1k, Qt, gamma, alpha, xd, Qm);
                end
                k = k+1; %increment the simulation step
                if mod(k, P) == 0 && k > Ninit
                    Qt.trainedNet = SeriesNetwork(Qm.layers); % copy weights
                    J = cost_to_go(Qt, s);
                    dJnorm = [dJnorm norm(J - Jlast)];
                    kJ = [kJ, k];
                    Jlast = J;
                    vi_plot(J, [], theta_bins, thetadot_bins, fh);
                    figure(11)
                    plot(kJ, dJnorm);
                    title('Change in cost-to-go');%%
                    xlabel('Simulation step');
                    ylabel('Norm difference');
                    drawnow;
                end
            end
        end
    else
        %load the target network from data
        Qtt = load('QPendNet.mat');
        Qt = Qtt.Qt;
    end
    
    J = cost_to_go(Qt,s);% get the cost-to-go

    vi_plot(J, [], theta_bins, thetadot_bins, fh);
    
    x0 = [0; 0; 0; 0];
    epsilon = 0;
    xsave = zeros(4, 101);
    for k = 1:101
        draw_cartpole(0, x0, param);
        xsave(:,k) = x0;
        a_ind = get_action(Qt, x0, epsilon);
        a0k(:,k) = a_ind;
        x1 = x0 + cartpole_dynamics(x0, Qt.a(a_ind), param) * dt;
        x1 = normalize(x1, x_bins, theta_bins, xdot_bins, thetadot_bins);
        x0 = x1;
    end
        figure(1)
    subplot(2,1,1)
    plot(xsave(2,:));
    title('Theta over time');
    xlabel('Time step'); ylabel('Theta (rad)');
    subplot(2,1,2)
    plot(xsave(4,:));
    title('Angular velocity over time');
    xlabel('Time step'); ylabel('Theta dot (rad/s)');
end

% ===================================================================
% Train the main network
% ===================================================================
function Qm = train(x0k, a0k, x1k, Qt, gamma, alpha, xd, Qm)
    inds = randi(length(a0k), 500, 1);
    xs = x0k(:, inds);
    as = a0k(:, inds);
    x1s = x1k(:, inds);
    % Compute TD targets using the target network.
    Qj = td_update(xs, as, x1s, Qt, gamma, alpha, xd);
    % Train the main network using these targets.
    Qm = learn_q_factor(xs, Qj, Qm);
end

% ===================================================================
% Update the Q-values using the TD error.
% ===================================================================
function Qxnew = td_update(x0, a_ind, x1, Qt, gamma, alpha, xd)
    Qx = predict(Qt.trainedNet, x0');  % [batch x N_actions]
    Qx1 = predict(Qt.trainedNet, x1');
    Qhat0 = get_q_factor(a_ind, Qx);
    Qhat1 = get_Qstar(Qx1);
    u = Qt.a(a_ind);
    c = cost_function(x0, u, xd);
    c = c(:);
    Qnew = (1 - alpha) * Qhat0 + alpha * (c + gamma * Qhat1);
    Qxnew = set_q_factor(a_ind, Qnew, Qx);
end

function Q = learn_q_factor(xs, Qnew, Q)
    trainedNet = trainNetwork(xs', Qnew, Q.layers, Q.opts);
    Q.trainedNet = trainedNet;
    Q.layers = trainedNet.Layers;
end

% ===================================================================
% Compute the cost-to-go from the 
% ===================================================================
function J = cost_to_go(Q, s)
    Qx = predict(Q.trainedNet, s');
    [J, ~] = min(Qx, [], 2);
end

% ===================================================================
% Dynamics model for the pendulum 
% ===================================================================
function xdot = cartpole_dynamics(x, u, param)
    mc = param.mc;
    mp = param.mp;
    l = param.l;
    g = param.g;
    b = param.b;
    d = param.d;
    x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
    s = sin(x2); c = cos(x2);
    xdot = zeros(4,1);
    xdot(1) = x3;
    xdot(2) = x4;
    xdot(3) = (u - b*x3 + d*x4*c/l + mp*s*(l*x4^2 + g*c)) / (mc + mp*s^2);
    xdot(4) = (-u*c + b*x3*c - d*(mc+mp)*x4/(mp*l) - mp*l*x4*c*s - (mc+mp)*g*s/(mp*l)) / (l*(mc+mp*s^2));
end

% ===================================================================
% This function defines the instantaneous cost (i.e. g(x,u))
% Note that X and u are vectors
% ===================================================================
function C = cost_function(X, u, xd)
    global dt;
    [Q, R] = get_QR;
    e = X - xd;
    C = Q(1,1)*(e(1,:).^2) + Q(2,2)*(e(2,:).^2) + ...
        Q(3,3)*(e(3,:).^2) + Q(4,4)*(e(4,:).^2) + R*u.^2;
end

function [Q, R] = get_QR
    global dt;
    Q = diag([1, 500, 20, 50]) * dt;
    R = 1 * dt;
end

% ===================================================================
% Compute the policy given from a particular state using the optimal
% q-factor (IMPLEMENT EPSILON-GREEDY APPROACH)
% ===================================================================
function ind = get_action(Q, x, epsilon)
    if rand < epsilon
        ind = randi(length(Q.a));
    else
        Qvalues = predict(Q.trainedNet, x');
        [~, ind] = min(Qvalues);
    end
end

% ===================================================================
% Compute the optimal Q-value from all actions.
% ===================================================================
function Qstar = get_Qstar(Qx)
    Qstar = min(Qx, [], 2);
end

% ===================================================================
% Get the Q-factor for set of actions (FILL IN)
% ===================================================================
function Q0 = get_q_factor(as, Qx)
    N = size(Qx,1);
    idx = sub2ind(size(Qx), (1:N)', as(:));
    Q0 = Qx(idx);
end

% ===================================================================
% Set the q-factor at a particular state given a set of actions (FILL IN)
% ===================================================================
function Qx = set_q_factor(as, Q0, Qx)
    N = size(Qx, 1);
    idx = sub2ind(size(Qx), (1:N)', as(:));
    Qx(idx) = Q0;
end

% ==============================================================
% This is the draw function for the brick
% ==============================================================
function draw_cartpole(t, x, param)
    l = param.l;
    persistent hFig base a1 raarm wb lwheel;
    % Reinitialize hFig if it's empty or not a valid figure handle
    if isempty(hFig) || ~ishandle(hFig)
        hFig = figure(25);
        set(hFig,'DoubleBuffer','on');
        a1 = l + 0.25;
        av = pi*[0:.05:1];
        theta_vals = pi*[0:0.05:2];
        wb = 0.3; hb = 0.15;
        aw = 0.01;
        wheelr = 0.05;
        angles = 0:0.05:2*pi;
        lwheel = [-wb/2 + wheelr*cos(angles); -hb - wheelr + wheelr*sin(angles)]';
        base = [wb*[1 -1 -1 1]; hb*[1 1 -1 -1]]';
        arm = [aw*cos(av-pi/2) -a1 + aw*cos(av+pi/2)
               aw*sin(av-pi/2) aw*sin(av+pi/2)]';
        raarm = [(arm(:,1).^2+arm(:,2).^2).^.5, atan2(arm(:,2),arm(:,1))];
    end

    figure(hFig); cla; hold on; view(0,90);
    % Draw cart as a rectangle.
    patch(x(1) + base(:,1), base(:,2), 0*base(:,1), 'b', 'FaceColor', [0.3 0.6 0.4]);
    % Draw wheels.
    patch(x(1) + lwheel(:,1), lwheel(:,2), 0*lwheel(:,1), 'k');
    patch(x(1) + wb + lwheel(:,1), lwheel(:,2), 0*lwheel(:,1), 'k');
    % Draw pole.
    patch(x(1) + raarm(:,1).*sin(raarm(:,2) + x(2) - pi), ...
          -raarm(:,1).*cos(raarm(:,2) + x(2) - pi), ...
          1+0*raarm(:,1), 'r', 'FaceColor', [0.9 0.1 0]);
    % Draw pivot.
    plot3(x(1) + l*sin(x(2)), -l*cos(x(2)), 1, 'ko', 'MarkerSize',10,'MarkerFaceColor','b');
    plot3(x(1), 0, 1.5, 'k.');
    title(['t = ', num2str(t, '%.2f'), ' sec']);
    set(gca, 'XTick',[],'YTick',[]);
    axis image; axis([-2.5 2.5 -2.5*l 2.5*l]);
    drawnow;
end

% ===================================================================
% Plot the cost-to-go over the (theta, thetadot) grid.
% ===================================================================
function vi_plot(J, ~, theta_bins, thetadot_bins, fh)
    set(0, 'CurrentFigure', fh);
    clf;
    n1 = length(theta_bins); n2 = length(thetadot_bins);
    imagesc(theta_bins, thetadot_bins, reshape(J, n1, n2)');
    axis xy; colorbar;
    title('Cost-to-go (Theta vs Theta dot)');
    xlabel('Theta (rad)'); ylabel('Theta dot (rad/s)');
    drawnow;
end

% ===================================================================
% ===================================================================
function s = normalize(s, x_bins, theta_bins, xdot_bins, thetadot_bins)
    N = size(s,2);
    % Wrap theta (state index 2)
    s(2,:) = mod(s(2,:), 2*pi);
    % Clamp cart position (state index 1)
    s(1,:) = min(max(s(1,:), min(x_bins)), max(x_bins));
    % Clamp cart velocity (state index 3)
    s(3,:) = min(max(s(3,:), min(xdot_bins)), max(xdot_bins));
    % Clamp pole angular velocity (state index 4)
    s(4,:) = min(max(s(4,:), min(thetadot_bins)), max(thetadot_bins));
end
