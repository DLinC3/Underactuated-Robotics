function grid_world_sarsa_qlearning
    close all
    clear all

    %grid world
    x=1:1:8;
    y=1:1:4;
    
    [X, Y] =meshgrid(x,y);
    Z = 0*X;

    use_pits = true;

    %add the obstacle positions
    xyobs{1} = [2,2];
    xyobs{2} = [3,1];
    xyobs{3} = [3,2];

    if use_pits==true
        xypits{1} = [6,3];
        xypits{2} = [6,4];    
    else
        xyobs{4} = [6,2];
        xyobs{5} = [6,3];
        xyobs{6} = [6,4];
        xyobs{7} = [7,2];
        xypits=[];
    end

    xd = [8,4];

    %create the grid world
    grid_world = create_grid_world(X,Y,Z,xyobs,xypits);
    %plot_grid_world(grid_world)

    [m n] = size(Z);
    Q=0*repmat(randn(m,n),1,1,5);

    Qlast = Q;

    Qc.X=X;%q-function x-grid points
    Qc.Y=Y;%q-function y-grid points
    Qc.Q=Q;%q-function

    [m n] = size(Z);

    %desired final state
    grid_world.xd = xd;

    %discout factor
    gamma = 0.9;

    N = 40;% number of episode time steps
    alpha = 0.5; %learning rate
    % 0 SARSA 1 Q-learning
    use_qlearning = 1;
    i=1;
    M = 500;%number of iterations (fill this in)
    J_norms = zeros(M, 1)
    for j=1:M
        j
        x0 = getRandomState(Qc, grid_world, xd); %get random initial state (fill this in)
        % epsilon = 0.2;%epsilon for eps-greedy approach
        epsilon = 1/(j+2);
        for k=1:N
                xsave(:,k) = x0;
                a0= get_action(Qc,x0,epsilon); %get the action
                x1=dynamics(x0,a0,grid_world);  %simulate forward
                a1= get_action(Qc,x1,epsilon); %get the action at the next time-step
                Qc = td_update(x0,a0,x1,a1,Qc,gamma,alpha,xd,grid_world,use_qlearning);%temporal difference update function
                x0 = x1;
        end
        J = cost_to_go(Qc);% get the cost-to-go
        J_norms(j) = norm(J(:));
        figure(1)
        hold on
        imagesc(X(1,:), Y(:,1), J)
        plot(xsave(1,:),xsave(2,:),'r')% plot the cost-to-go
        colorbar
        grid on
        drawnow;
        hold off      

        Q = Qc.Q;
        eQ = abs(Q(:)-Qlast(:));%look at the q-function error between iterations
        normQ(i) = norm(Q(:));
        i = i + 1;
        eQ = max(eQ);
        Qlast = Qc.Q;
      
    end
    figure(1)
    hold on
    plot_obstacles(grid_world) %plot the ostacles
    hold off

%     figure(1)
%     J= J(:);
%     for i=1:length(J)
%         text(X(i),Y(i),num2str(J(i),2))
%     end

    figure(2)
    plot(normQ); %plot the iteration Q-function norm error

    x0=[1;1];%intial conditions
    epsilon = 0;
    %simulate the policy from the initial conditions
   for k=1:100
        xsave(:,k) = x0;
        a0= get_action(Qc,x0,epsilon);
        x1=dynamics(x0,a0,grid_world);  
        x0 = x1;
   end
   figure(1)
   hold on
   plot(xsave(1,:),xsave(2,:),'*-g')%plot the trajectory for the optimal policy
   hold off
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
function grid_world = create_grid_world(X,Y,Z,xyobs,xypits)
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
    P=Z*0;
    for i=1:length(xypits)
        P(find_index(xypits{i}(1),xypits{i}(2),X,Y)) = 1;
    end
    grid_world.P = P;
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

function ret = check_pits(x,grid_world)
    inds = find(grid_world.P>0);
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
% Update the Q-function using the temporal difference error (fill this in)
% ===================================================================
function Qc = td_update(x0,a0,x1,a1,Qc,gamma,alpha,xd, grid_world, use_qlearning)
    % 获取当前Q值
    current_Q = get_q_factor(x0, a0, Qc);
    
    % 计算即时成本（处理目标状态和陷阱）
    g = cost_function(x0, xd, grid_world);
    
    % 判断是否为终止状态
    is_terminal = isequal(x0,xd) || check_pits(x0,grid_world);
    
    if is_terminal
        target = g;
    else
        if use_qlearning
            % Q-learning: 选择下一状态的最小Q值（因是最小化问题）
            q_values = arrayfun(@(a) get_q_factor(x1,a,Qc), 0:4);
            target = g + gamma * min(q_values);
        else
            % SARSA: 使用实际采取的下一动作Q值
            next_Q = get_q_factor(x1,a1,Qc);
            target = g + gamma * next_Q;
        end
    end
    
    % 执行Q值更新
    updated_Q = current_Q + alpha * (target - current_Q);
    Qc = set_q_factor(x0, a0, updated_Q, Qc);
end

% ===================================================================
% Get a random state
% ===================================================================
function x = getRandomState(Qc,grid_world,xd)
    ret = 1;
    x = xd;
    while (ret==1) || (norm(x-xd)==0)
        randidx = randi(length(Qc.X(:)));
        x(1) = Qc.X(randidx);
        x(2) = Qc.Y(randidx);
        ret=check_occupied(x,grid_world);
    end
end


% ===================================================================
% Get the Q-factor for state x and action a
% ===================================================================
function Q0= get_q_factor(x,a,Q)
    [i,j0]=find(Q.X==x(1));
    [i,j]=find(Q.Y(:,j0(1))==x(2));
    Q0=Q.Q(i,j0(1),a+1);
end

% ===================================================================
% Set the q-factor for state x and action a
% ===================================================================
function Q= set_q_factor(x,a,Q0,Q)
    [i,j0]=find(Q.X==x(1));
    [i,j]=find(Q.Y(:,j0(1))==x(2));
    Q.Q(i,j0(1),a+1)=Q0;
end

% ===================================================================
% compute the cost-to-go from the q-factor
% ===================================================================
function J = cost_to_go(Q)

    [J, ind]= min(Q.Q,[],3);
   
end

% ===================================================================
% Define the one step cost.
% ===================================================================
function g = cost_function(x,xd,grid_world)
    if(x(1)==xd(1) && x(2)==xd(2))
        g=0;
    elseif check_pits(x,grid_world)==1  
        g=10;
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
        elseif check_pits(x,grid_world) == 1
            x = x;
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
% q-factor (fill this in)
% ===================================================================
function a= get_action(Q,x,epsilon)
    % 获取所有动作的Q值
    q_values = arrayfun(@(a) get_q_factor(x,a,Q), 0:4);
    
    % 以epsilon概率随机探索
    if rand < epsilon
        a = randi(5)-1; % 随机选择0-4
    else
        % 选择最小Q值动作（因是最小化问题）
        [~, idx] = min(q_values);
        a = idx-1; % 转换为0-based索引
        
        % 处理多个最优动作的情况
        min_val = min(q_values);
        best_actions = find(q_values == min_val);
        a = best_actions(randi(length(best_actions))) - 1;
    end
end