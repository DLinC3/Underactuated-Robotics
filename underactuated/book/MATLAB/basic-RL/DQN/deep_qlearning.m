function deep_qlearning
    global dt;
    dt = 0.1;
    xd = [pi,0];
    
    %discout factor
    gamma = 0.8;

    %define the discrete actions
    Qm. a = linspace(-2,2,10);
    
    %% Define a layer architecture for Q function (FILL IN)
    layers = [ ...
        featureInputLayer(2, "Name", "myFeatureInputLayer")
        fullyConnectedLayer(64, "Name", "myFullyConnectedLayer1")
        reluLayer("Name", "myReLu1")
        fullyConnectedLayer(32, "Name", "myFullyConnectedLayer2")
        reluLayer("Name", "myReLu2")
        fullyConnectedLayer(length(Qm.a), "Name", "myFullyConnectedLayer3")%10
        regressionLayer("Name", "myRegressionLayer")
    ];
    %% Define options for the training
    opts = trainingOptions('adam', ...
        'MaxEpochs',100, ...
        'InitialLearnRate',0.01,...
        'MiniBatchSize',1024, ...
        'Verbose',false, ...
        Plots="none");
    % create two networks, one for training and one for evaluation (target)
    Qm.layers = layers; % main network for training (main network)
    Qm.opts = opts;
    Qt =Qm; %network for evaluation (target network)

    % create a mesh for visualizing the cost-to-go    
    q_bins = linspace(0,2*pi,25);   
    qdot_bins = linspace(-10,10,25);
    [q qdot] = ndgrid(q_bins,qdot_bins);
    s = [reshape(q,1,numel(q)); reshape(qdot,1,numel(qdot))];    
    
    N = 20; %train the network every "N" simulation steps
    Ninit = 1000; %number of simulation steps before training starts
    Nepisode = 200; % maximum number of steps in an episode;
    %initialize the networks
    x0k = zeros(2,N);
    a0k = zeros(1,N);
    Qmk = zeros(N,length(Qm.a));
    Qm= learn_q_factor(x0k,Qmk,Qm);
    Qt= learn_q_factor(x0k,Qmk,Qm);
    fh = figure(10);
    Jlast = zeros(length(s),1); %cost function from last iteration
    loadFromData =false;
    dJnorm = [];
    kJ = [];
    
    if loadFromData==false
         
        alpha = 0.8; %learning rate
        M = 300;%number of episodes
        P = 100; %number of steps before copying weights to target network
        k = 1;
        for i = 1:M

            x0=[0;0];
            
            done = false;
            epsilon = 0.01 + (1 - 0.01) * exp(-0.01 * i); %update epsilon with every episode
            j = 1;
            % while the episode is not complete
            while done==false 
                display('epsiode || steps in episode || simulation steps')
                [i j k]
                draw(0,x0) %draw the robot using the current state
                x0k(:,k) = x0; %add the state to the history of data
                a_ind = get_action(Qt,x0,epsilon); %get the action index using epsion-greedy approach
                a0k(:,k) = a_ind; %save the action index
                x1=x0+dynamics(x0,Qt.a(a_ind))*dt;  %simulate the dynamics forward
                if (x1(2)>max(qdot_bins)) || (x1(2)<min(qdot_bins)) %if we exceed the velocity bounds, the game is over
                    done = true;
                elseif(j>Nepisode) %if we exceed the episode time horizon, the game is over
                    done = true;
                else
                end
                j=j+1;
                % account for angle wrapping and velocity bounds
                x1 = normalize(x1,q_bins,qdot_bins);
                % save the next state
                x1k(:,k)= x1;
                % set 
                x0 = x1;
                
                %train the main network every N simulation steps
                if mod(k,N)==0 && k>Ninit
                    Qm = train(x0k,a0k,x1k,Qt,gamma,alpha,xd,Qm);
                end
                k = k+1; %increment the simulation step
                
                if mod(k,P)==0 && k>Ninit
                    Qt.trainedNet= SeriesNetwork(Qm.layers); %copy the weights from the main network to the target network
                    J = cost_to_go(Qt,s);% get the cost-to-go at the grid points
                    dJnorm = [dJnorm norm(J-Jlast)]; %compare the current cost-to-go with the prior cost-to-go
                    kJ =[kJ k]; % save the simulation step
                    Jlast = J;
                    vi_plot(J,x1k,q_bins,qdot_bins,fh);  %plot the cost-to-go  
                    figure(11)
                    plot(kJ,dJnorm); %plot the change in the cost-to-go
                end
            
            end
        
        end
    else
        %load the target network from data
        Qtt = load('QPendNet.mat');
        Qt = Qtt.Qt;
    end
    
    J = cost_to_go(Qt,s);% get the cost-to-go
    vi_plot(J,[],q_bins,qdot_bins,fh);  
    
    
    x0=[0;0];%intial conditions
    epsilon = 0;
    %simulate the policy from the initial conditions
    for k=1:100
        draw(0,x0)
        xsave(:,k) = x0;
        a_ind= get_action(Qt,x0,epsilon); %get the action
        a0k(:,k) = a_ind;
        x1=x0+dynamics(x0,Qt.a(a_ind))*dt; %simulate the dynamics forward
        x1 = normalize(x1,q_bins,qdot_bins);%apply angle wrapping
        x0 = x1;
    end
    figure(1)
    hold on
    plot(xsave(1,:),xsave(2,:),'*-g') %plot the trajectory for the optimal policy
    hold off
end
% ===================================================================
% Train the main network
% ===================================================================
function Qm = train(x0k, a0k, x1k, Qt, gamma, alpha, xd, Qm)
    inds = randi(length(a0k), 500, 1);
    xs = x0k(:, inds);
    as = a0k(:, inds);
    x1s = x1k(:, inds);
    Qj = td_update(xs, as, x1s, Qt, gamma, alpha, xd);
    Qm = learn_q_factor(xs, Qj, Qm);
end

% ===================================================================
% Update the Q-function using the temporal difference error (FILL IN)
% ===================================================================
function Qxnew = td_update(x0, a_ind, x1, Qt, gamma, alpha, xd)
    Qx = predict(Qt.trainedNet, x0');  % size: [batch_size x N_actions]
    Qx1 = predict(Qt.trainedNet, x1');
    
    Qhat0 = get_q_factor(a_ind, Qx);
    Qhat1 = get_Qstar(Qx1);
    
    u = Qt.a(a_ind);
    
    c = cost_function(x0, u, xd);
    c = c(:);  % Convert to column vector
    
    Qnew = (1 - alpha) * Qhat0 + alpha * (c + gamma * Qhat1);%td error
    
    Qxnew = set_q_factor(a_ind, Qnew, Qx);
end




function Q = learn_q_factor(xs, Qnew, Q)% train the Q function (FILL IN)
    trainedNet = trainNetwork(xs', Qnew, Q.layers, Q.opts);
    Q.trainedNet = trainedNet;
    Q.layers = trainedNet.Layers;
end



% ===================================================================
% compute the cost-to-go from the q-factor
% ===================================================================
function J = cost_to_go(Q,s)
    Qx = predict(Q.trainedNet,s');
    [J, ind]= min(Qx,[],2);
end

% ===================================================================
% Dynamics model for the pendulum 
% ===================================================================
function xdot = dynamics(x,u)
 m = 1;   % kg
    l = .5;  % m
    b = 0.1; % kg m^2 /s
    lc = .5; % m
    I = .25; %m*l^2; % kg*m^2
    g = 9.8;
    xdot = [x(2,:); (u - m*g*l*sin(x(1,:)) - b*x(2,:))/I];

end

% ===================================================================
% This function defines the instantaneous cost (i.e. g(x,u))
% Note that X and u are vectors
% ===================================================================
function C = cost_function(X,u,Xd)
    global dt
    [Q,R]=get_QR;
     
    C=(Q(1,1)*(X(1,:)-pi).^2 + Q(2,2)*X(2,:).^2)+R*u.^2;
end

function [Q,R] = get_QR
    global dt;
    
    Q = diag([10 1]).*dt; % <== the dt just discretizes
    R = 1*dt;

end

% ===================================================================
% Compute the policy given from a particular state using the optimal
% q-factor (IMPLEMENT EPSILON-GREEDY APPROACH)
% ===================================================================
function ind = get_action(Q, x, epsilon)
    if rand < epsilon
        % select a random action index
        ind = randi(length(Q.a));
    else
        Qvalues = predict(Q.trainedNet, x');  % output is a row vector of Q-values
        [~, ind] = min(Qvalues);
    end
end

% ===================================================================
% Compute the optimal q-factor
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
    function draw(t,x)
      % Draw the pendulum.  
      persistent hFig base a1 a2 ac1 ac2 raarm t0;

      if (isempty(hFig))
        hFig = figure(25);
        set(hFig,'DoubleBuffer', 'on');
        
        a1 = 0.75;  ac1 = 0.415;
        av = pi*[0:.05:1];
        rb = .03; hb=.07;
        aw = .01;
        base = rb*[1 cos(av) -1 1; -hb/rb sin(av) -hb/rb -hb/rb]';
        arm = [aw*cos(av-pi/2) -a1+aw*cos(av+pi/2)
          aw*sin(av-pi/2) aw*sin(av+pi/2)]';
        raarm = [(arm(:,1).^2+arm(:,2).^2).^.5, atan2(arm(:,2),arm(:,1))];
      end
            
      figure(hFig); cla; hold on; view(0,90);
      patch(base(:,1), base(:,2),1+0*base(:,1),'b','FaceColor',[.3 .6 .4])
      patch(raarm(:,1).*sin(raarm(:,2)+x(1)-pi),...
        -raarm(:,1).*cos(raarm(:,2)+x(1)-pi), ...
        0*raarm(:,1),'r','FaceColor',[.9 .1 0])
      plot3(ac1*sin(x(1)), -ac1*cos(x(1)),1, 'ko',...
        'MarkerSize',10,'MarkerFaceColor','b')
      plot3(0,0,1.5,'k.')
      title(['t = ', num2str(t(1),'%.2f') ' sec']);
      set(gca,'XTick',[],'YTick',[])
      
      axis image; axis([-1.0 1.0 -1.0 1.0]);
      drawnow;
    end  

% ===============================================================
% This function plots the value function and policy
%================================================================
function vi_plot(J,xk,q_bins,qdot_bins,fh)
set(0, 'CurrentFigure', fh);
clf reset;
n1 = size(q_bins,2); n2 = size(qdot_bins,2);
imagesc(q_bins,qdot_bins,reshape(J,n1,n2)'); axis xy;
%plot(xk(1,:),xk(2,:),'*k')
%hold off
drawnow; 
end

% ==============================================================
% This function imposes limits on the state
% ==============================================================
function s = normalize(s,q_bins,qdot_bins)
% impose limits
    N = size(s,2);
    s(1,:) = mod(s(1,:),2*pi);
    smax = repmat([q_bins(end);qdot_bins(end)],1,N); ind = s>smax;
    s(ind) = smax(ind);
    smin = repmat([q_bins(1);qdot_bins(1)],1,N); ind = s<smin;
    s(ind) = smin(ind);
end