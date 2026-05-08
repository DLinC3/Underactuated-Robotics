% ===================================================================
% Function implements infinite time value iteration for grid world
% ===================================================================
function grid_world_value_iteration_stochastic
    close all
    clear all

    %create the 4x8 grid
    xx=1:1:8;
    yy=1:1:4;    
    [X, Y] =meshgrid(xx,yy);
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
    Jlast =J;
    Xi=X(:);
    Yi=Y(:);
    Ji=J(:);
    Jc.X=X;
    Jc.Y=Y;
    Jc.J=J;

    e=1e6;
    % figure(1)
    hold on
    j=0;

    %set the desired position
    xd = [8,4];
    %set the discount factor
    gamma = .9;

    %loop through until the cost-to-go has converged
    while(e>0.01)
    j=j+1
        for i=1:length(Ji)
            % check if state is occupied    
            ret = check_occupied([Xi(i); Yi(i)],grid_world);
            if ret==0
            Ji(i) = bellman_backup([Xi(i); Yi(i)],Jc,grid_world,xd,gamma);
            else
            
            end
        end
    [m,n]= size(Z);    
    J=reshape(Ji,m,n)
    e=norm(J(:)-Jlast(:));% error between the cost-to-go from one iteration to the next
    normJ(j) = norm(J(:));
    Jlast = J;
    Jc.J=J;

    %plot the cost-to-go
    figure(1)
    title("cost-to-go")
    hold on
    imagesc(X(1,:), Y(:,1), J)
    colorbar
    grid on
    drawnow;
    hold off
    end
    
    %plot obstacles
    figure(1)
    hold on
    plot_obstacles(grid_world)
    hold off

    x=[1 4];
    %compute the policy from state s
    [pi, xs]= get_policy(x,Jc,grid_world,xd,gamma);
    figure(1)
    hold on
    plot(xs(:,1),xs(:,2),'r')
    hold off

    figure(1)
    J= J(:);
    for i=1:length(J)
        text(X(i),Y(i),num2str(J(i),2))
    end

    figure(2)
    plot(normJ(2:end))
    title("Convergence of norm(J)")
    xlabel("iteration")
    ylabel("norm(J)")
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

% ===================================================================
% plot grid world
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
% This function applies the bellman back-up from state x and finds the 
% minimum cost-to-go from state x given the cost-to-go at x'
% ===================================================================
function J0 = bellman_backup(x, J, grid_world, xd, gamma)
min_total = inf;
% Iterate over all possible actions (0-4)
for a = 0:4
    expected_J = 0;
    % Calculate expected cost for each possible executed action
    for a_execute = 0:4
        % Determine transition probability
        if a_execute == a
            p = 0.6;
        else
            p = 0.1;
        end
        xi = dynamics_nominal(x, a_execute, grid_world);
        J_xi = cost_to_go(xi, J);
        expected_J = expected_J + p * J_xi;
    end
    g = cost_function(x, xd);
    total_cost = g + gamma * expected_J;
    if total_cost < min_total
        min_total = total_cost;
    end
end
J0 = min_total;
end
% ===================================================================
% Get the cost-to-go from state x
% ===================================================================
function J0= cost_to_go(x,J)
    [i,j0]=find(J.X==x(1));
    [i,j]=find(J.Y(:,j0(1))==x(2));
    J0=J.J(i,j0(1));
end

% ===================================================================
% Define the one step cost.
% ===================================================================
function g = cost_function(x,xd)
    if(x(1)==xd(1) && x(2)==xd(2))
        g=-1;
    else
        g=1;
    end
end

% ===================================================================
% dynamics_nominal model for grid world. 
% ===================================================================
function x = dynamics_nominal(x,u,grid_world)
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

% ===================================================================
% Keep state in grid world bounds
% ===================================================================
function x = applyBounds(x,occ_grid)
    if(x(1)<occ_grid.min_x)
        x(1)=occ_grid.min_x;
    end
    if(x(1)>occ_grid.max_x)
        x(1)=occ_grid.max_x;
    end
    if(x(2)<occ_grid.min_y)
        x(2)=occ_grid.min_y;
    end
    if(x(2)>occ_grid.max_y)
        x(2)=occ_grid.max_y;
    end
end

% ===================================================================
% Compute the policy given from a particular state using the optimal
% cost-to-go
% ===================================================================
function [pi, xs] = get_policy(x, J, grid_world, xd, gamma)
pi = [];
xs = [x];
while true
    Js = zeros(1,5);
    for a = 0:4
        expected_J = 0;
        for a_execute = 0:4
            if a_execute == a
                p = 0.6;
            else
                p = 0.1;
            end
            xi = dynamics_nominal(x, a_execute, grid_world);
            J_xi = cost_to_go(xi, J);
            expected_J = expected_J + p * J_xi;
        end
        g = cost_function(x, xd);
        Js(a+1) = g + gamma * expected_J;
    end
    [~, ind] = min(Js);
    best_action = ind-1;
    xi = dynamics_nominal(x, best_action, grid_world);
    x = xi;
    pi = [pi, best_action];
    xs = [xs; x];
    if all(x == xd)
        break;
    end
end
end