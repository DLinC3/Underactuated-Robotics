function grid_world_td_learning
    close all
    clear all

    x=1:1:8;
    y=1:1:4;
    
    [X, Y] =meshgrid(x,y);
    Z = 0*X;

    %add the obstacle positions
    xyobs{1} = [2,2];
    xyobs{2} = [3,1];
    xyobs{3} = [3,2];
    xyobs{4} = [6,2];
    xyobs{5} = [6,3];
    xyobs{6} = [6,4];
    xyobs{7} = [7,2];

    %create the grid world
    grid_world = create_grid_world(X,Y,Z,xyobs);
    %plot_grid_world(grid_world)

    J=0*Z;
    Jlast = J;

    Xi=X(:);
    Yi=Y(:);
    Ji=J(:);

    %the cost-to-go and its grid points
    Jc.X=X;
    Jc.Y=Y;
    Jc.J=J;

    [m n] = size(Z);

    %desired final state
    xd = [8,4];
    grid_world.xd = xd;

    %discout factor
    gamma = 0.9;

    figure(1)
    hold on
    N = 20; %episode horizon
    M=500;% number of iterations (fill this in)
    alpha = 0.1;%learnng rate (fill this in)
    for j=1:M
        j
        x0 = getRandomState(Jc, grid_world, xd); % fill in initial state
        for k=1:N
            k;
            xsave(:,k) = x0;
            a0= get_action(Jc,grid_world,x0,xd,gamma); %get the action for the current state
            x1=dynamics(x0,a0,grid_world);  %forwar simulate the dynamics
            Jc = td_update(x0,a0,x1,Jc,gamma,alpha,xd); %use the td-error to update the cost-to-go
            x0 = x1;
        end

        figure(2)
        hold on
        imagesc(X(1,:), Y(:,1), Jc.J) %plot the cost-to-go
        plot(xsave(1,:),xsave(2,:),'r') %plot the state trajectory
        colorbar
        grid on
        drawnow;
        hold off      

        J = Jc.J;
        J_norms(j) = norm(J(:)); % add J norms
        eJ = abs(J(:)-Jlast(:));
        eJ = max(eJ);
        Jlast = Jc.J;
      
    end
    figure(2)
    hold on
    plot_obstacles(grid_world)
    hold off

    figure(2)
    J= J(:);
    for i=1:length(J)
        text(X(i),Y(i),num2str(J(i),2))
    end

    x0=[1;1];%intial conditions
    %simulate the policy from the initial conditions
   for k=1:100
        xsave(:,k) = x0;
        a0= get_action(Jc,grid_world,x0,xd,gamma);
        x1=dynamics(x0,a0,grid_world);  
        x0 = x1;
   end
   figure(2)
   hold on
   plot(xsave(1,:),xsave(2,:),'*-g')%plot the trajectory for the optimal policy
   hold off

    figure(3);
    plot(1:M, J_norms);
    xlabel('Iteration');
    ylabel('2-norm of Cost-to-go');
    title('Convergence of Cost-to-go');
    grid on;
end

% ===================================================================
% Plot the grid world
% ===================================================================
function plot_grid_world(grid_world)
    xs = grid_world.X(1,:);
    ys = grid_world.Y(:,1);
    dx = grid_world.dx;
    dy = grid_world.dy;
    xg = linspace(xs(1)-dx,xs(end)+dx,length(xs)+1);
    yg = linspace(ys(1)-dy,ys(end)+dy,length(ys)+1);
    hi = imagesc(xs,ys,grid_world.Z); 
    set(gca,'YDir','normal');
    hold on
    hm = mesh(xg,yg,zeros([length(ys) length(xs)]+1));
    hm.FaceColor = 'none';
    hm.EdgeColor = 'k';
end

% ===================================================================
% plot obstacles
% ===================================================================
function plot_obstacles(grid_world)
     inds = find(grid_world.Z>0);
     dx = 0.5;
     dy = 0.5;
     occxs=grid_world.X(inds);
     occys=grid_world.Y(inds);
     for i=1:length(inds)
        xp=[occxs(i)-dx occxs(i)-dx occxs(i)+dx occxs(i)+dx];
        yp=[occys(i)-dy occys(i)+dy occys(i)+dy occys(i)-dy];
        patch(xp,yp,[1 1 1]*.3,'LineStyle','None')
     end
end

function k = find_index(x,y,X,Y)
    inds = find(X==x(1));
    ind = find(Y(inds)==y);
    k = inds(ind);
end

% ===================================================================
% Create the grid world
% ===================================================================
function grid_world = create_grid_world(X,Y,Z,xyobs)
    grid_world.X = X;
    grid_world.Y = Y;
    grid_world.Z = Z;
    xs = grid_world.X(1,:);
    ys = grid_world.Y(:,2);
    grid_world.dx = abs(xs(2)-xs(1))/2;
    grid_world.dy = abs(ys(2)-ys(1))/2; 
    grid_world.max_x = max(X(:));
    grid_world.max_y = max(Y(:));
    grid_world.min_x = min(X(:));
    grid_world.min_y = min(Y(:));
    for i=1:length(xyobs)
        Z(find_index(xyobs{i}(1),xyobs{i}(2),X,Y)) = 1;
    end
    grid_world.Z = Z;
end

% ===================================================================
% Check if state is occupied
% ===================================================================
function ret = check_occupied(x,grid_world)
    inds = find(grid_world.Z>0);
    occxs=grid_world.X(inds);
    occys=grid_world.Y(inds);
    ret = 0;
    for i=1:length(occxs)
        if (x(1)==occxs(i))
                if (x(2)==occys(i))
                    ret = 1;
                end
        end
    end
end

% ===================================================================
% The temporal difference update for the cost-to-go (fill this in)
% ===================================================================
function Jc = td_update(x0,a0,x1,Jc,gamma,alpha,xd)
    % Get cost estimates for the current state and next state
    J0 = cost_to_go(x0, Jc);
    J1 = cost_to_go(x1, Jc);
    
    % Calculate the immediate cost 
    % (the target state cost is 0, the others are 1)
    if isequal(x0(:),xd(:))
        g = 0;
    else
        g = 1;
    end
    
    % TD update
    td_target = g + gamma * J1;
    td_error = td_target - J0;
    
    % update
    J0_updated = J0 + alpha * td_error;
    Jc = set_cost_to_go(x0, J0_updated, Jc);
end

% ===================================================================
% Get a random initial state
% ===================================================================
function x = getRandomState(Jc,grid_world,xd)
    ret = 1;
    x=xd;
    while (ret==1) || (norm(x-xd)==0)
        randidx = randi(length(Jc.J(:)));
        x(1) = Jc.X(randidx);
        x(2) = Jc.Y(randidx);
        ret=check_occupied(x,grid_world);
    end
end

% ===================================================================
% Get the cost-to-go from state x
% ===================================================================
function J0= cost_to_go(x,J)
    [i,j0]=find(J.X==x(1));
    [i,j]=find(J.Y(:,j0(1))==x(2));
    J0=J.J(i,j0(1));
end

function J= set_cost_to_go(x,J0,J)
    [i,j0]=find(J.X==x(1));
    [i,j]=find(J.Y(:,j0(1))==x(2));
    J.J(i,j0(1))=J0;
end


% ===================================================================
% Define the one step cost.
% ===================================================================
function g = cost_function(x,xd)
x
    if(x(1)==xd(1) && x(2)==xd(2))
        g=0;
    else
        g=1;
    end
end

% ===================================================================
% Dynamics model for grid world. 
% ===================================================================
function x = dynamics(x,u,grid_world)
    if x==grid_world.xd
        x = grid_world.xd;
    else
        switch(u)
            case 0
                x=x;        
            case 1
                x(1)=x(1)+1;
                ret = check_occupied(x,grid_world);
                if ret ==1
                    x(1) = x(1)-1;
                end
            case 2
                x(1)=x(1)-1;
                ret = check_occupied(x,grid_world);
                if ret ==1
                    x(1) = x(1)+1;
                end
            case 3
                x(2)=x(2)+1;
                ret = check_occupied(x,grid_world);
                if ret ==1
                    x(2) = x(2)-1;
                end
            case 4            
                x(2)=x(2)-1;
                ret = check_occupied(x,grid_world);
                if ret ==1
                    x(2) = x(2)+1;
                end
        end
    x = applyBounds(x,grid_world);
    end
end

% ===================================================================
% Keep state in grid world bounds
% ===================================================================
function x = applyBounds(x,grid_world)
    if(x(1)<grid_world.min_x)
        x(1)=grid_world.min_x;
    end
    if(x(1)>grid_world.max_x)
        x(1)=grid_world.max_x;
    end
    if(x(2)<grid_world.min_y)
        x(2)=grid_world.min_y;
    end
    if(x(2)>grid_world.max_y)
        x(2)=grid_world.max_y;
    end
end

% ===================================================================
% Compute the policy given from a particular state using the optimal
% cost-to-go (fill this in with the greedy policy)
% ===================================================================
function a= get_action(J,grid_world,x,xd,gamma)
    min_value = inf;
    best_actions = [];
    
    % Define all possible actions
    % (0:stay 1:right 2:left 3:up 4:down)
    actions = 0:4;
    
    % all actions to calculate Q value
    for action = actions
        % next state
        x_next = dynamics(x, action, grid_world);
        
        % calculate Q
        if isequal(x(:),xd(:))
            q_value = 0;
        else
            q_value = 1 + gamma * cost_to_go(x_next, J);
        end
        
        % optimal action
        if q_value < min_value
            min_value = q_value;
            best_actions = action;
        elseif q_value == min_value
            best_actions = [best_actions, action];
        end
    end
    
    
    a = best_actions(randi(length(best_actions)));
end