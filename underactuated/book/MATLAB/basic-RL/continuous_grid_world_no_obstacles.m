function continuous_grid_world_no_obstacles
    close all
    clear all
    
    %set up and create the grid world
    x=[1:.5:8];
    y=[1:.5:4];
    [X, Y] =meshgrid(x,y);
    Z = 0*X;
    xyobs = [];
    grid_world = create_grid_world(X,Y,Z,xyobs);
    
    %plot the grid world
    plot_grid_world(grid_world)
    
    n = length(Z(:));
    
    %set up the cost-to-go
    Ji=zeros(n,1);
    J.X=X;
    J.Y=Y;
    J.J=Z;
    J.type='linear';
    
    switch J.type
        case 'linear'
            J.thetas=[0;0;0];
        case 'quadratic'
            J.thetas=[0;0;0;0;0;0];
    end
    
    %variable to store the cost-to-go parameters
    thetas_last = J.thetas;
    
    x_base = 1:0.1:8;   
    y_base = 1:0.1:4;
    [X_base,Y_base] = meshgrid(x_base,y_base);

    Xi = [X_base(:)];
    Yi = [Y_base(:)];
    
    n = length(Xi);
    % ========================================================================

    
    %plot the sample points on the grid world
    % figure(5)
    % plot(Xi,Yi,'*r')
    % hold off
    
    maxIterations = 2000;%maximum number of iterations
    
    for j=1:maxIterations
        %loop through the sample points
        for i=1:n
            Ji(i) = bellman_backup([Xi(i); Yi(i)],J,grid_world);
        end

        %%% code for approximating the cost-to-go here%%%%%%%%%%%%%%%%%%%%%
        J = fit_cost_to_go(Xi, Yi, Ji, J);
        %check parameters for convergence
        J.thetas
        e=norm(J.thetas-thetas_last)
        if e<0.001
            break;
        end
        thetas_last = J.thetas;
    end
    
    %loop through the grid points and evaluate the cost-to-go
    for i=1:length(X(:))
        Jhat(i)= cost_to_go([X(i),Y(i)]',J);
    end
    [m,n]= size(Z);    
    Jhat=reshape(Jhat,m,n);
    
    %plot the cost to go evaluated at the grid points
    figure(1)
    title("cost-to-go")
    hold on
    imagesc(X(1,:), Y(:,1), Jhat)
    colorbar
    grid on
    drawnow;
    hold off
    
    %plot the cost-to-go using MATLAB surf
    switch J.type
    case 'linear'
    plot_linear_cost_to_go(J)
    case 'quadratic'
    plot_quadratic_cost_to_go(J)
    end
    
    %select an initial condition
    x=[1 1]';
    
    N = 100;
    %for N time steps, execute the policy
    for k=1:N
    xs(:,k) = x;
    Jk(k) = cost_function(x);
    u = get_policy(x,J,grid_world);
    x = dynamics(x,u,grid_world);
    end

    %Jk=sum(Jk);
    
    %plot the resulting policy on the cost-to-go
    figure(1)
    hold on
    plot(xs(1,:),xs(2,:),'*r')
    hold off

end

%==========================================================================
%plot the quadratic cost-to-go
%==========================================================================
function plot_quadratic_cost_to_go(J)
    x=[1:.1:8];
    y=[1:.1:4];
    [X,Y] = meshgrid(x,y);
    Z=J.thetas(1)*X.^2+J.thetas(2)*Y.^2+J.thetas(3)*X+J.thetas(4)*Y+J.thetas(5)*X.*Y+J.thetas(6);
    figure(2)
    hold on
    surf(X,Y,Z);
    colorbar
    hold off
    figure(3);
    surf(X, Y, Z, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); 
    colormap(parula); 
    view(40, 30);
    xlabel('X');
    ylabel('Y');
    zlabel('Cost-to-Go');
    title('Quadratic Cost-to-Go (3d)');
    colorbar;
    grid on;
end

%==========================================================================
%plot the linear cost-to-go
%==========================================================================
function plot_linear_cost_to_go(J)
    x=[1:.1:8];
    y=[1:.1:4];
    [X,Y] = meshgrid(x,y);
    Z=J.thetas(1)*X+J.thetas(2)*Y+J.thetas(3);
    figure(2)
    hold on
    surf(X,Y,Z);
    colorbar
    hold off
    figure(3)
    surf(X, Y, Z, 'FaceAlpha', 0.8,'EdgeColor', 'none');
    colormap(parula);
    shading interp;
    view(40, 30);
    xlabel('X');
    ylabel('Y');
    zlabel('Cost-to-Go');
    title('Linear Cost-to-Go (3d)');
    colorbar;
end

%==========================================================================
%plot grid world
%==========================================================================
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

%==========================================================================
%create grid world
%==========================================================================
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

%==========================================================================
% apply the bellman backup
%==========================================================================
function Ji=bellman_backup(x,J,grid_world)
    gamma = 0.9;
    U1=[1 -1 0 0];
    U2=[0 0 1 -1];
    for i=1:length(U1)
        u(1) = U1(i);
        u(2) = U2(i);
        xi=dynamics(x,u',grid_world);
        Jis(i) = gamma*cost_to_go(xi,J) + cost_function(x);
    end
    [Ji, ind]= min(Jis);
end

%==========================================================================
% evaluate the cost-to-go
%==========================================================================
function J0= cost_to_go(x,Jc)

switch Jc.type
    case 'linear'
    param=[];
    phi = linear(x,param);    
    case 'quadratic'
    param=[];
    phi = quadratic(x,param);    
end

J0 = Jc.thetas'*phi;
end

%==========================================================================
% function defines the stage cost 
%==========================================================================
function g = cost_function(x)
    xd =[8;4];
    dt = .1;
    if (abs(x(1)-xd(1))<=.5) && (abs(x(2)-xd(2))<=.5)
        g = 0;    
    else
        g = 1*dt;
    end
    % change in c
    g=norm(x-xd)^2;
end

%==========================================================================
% check if state is occupied
%==========================================================================
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

%==========================================================================
% simuate the dynamics forward by one time-step
%==========================================================================
function x = dynamics(x,u,occ_grid)
    dt = .1;
    x(1,:) = x(1,:) + u(1,:)*dt;
    x(2,:) = x(2,:) + u(2,:)*dt;
    
    ret = check_occupied(x,occ_grid);
    if ret == 1
        x(1,:) = x(1,:) - u(1,:)*dt;
        x(2,:) = x(2,:) - u(2,:)*dt;
    end
    
    x = applyBounds(x,occ_grid);    
end

%==========================================================================
% apply bounds to the grid world
%==========================================================================
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

%==========================================================================
% get the policy
%==========================================================================
function u = get_policy(x,J,grid_world)
    gamma = 0.9;
    U1=[1 -1 0 0];
    U2=[0 0 1 -1];
    for i=1:length(U1)
        u(1,i) = U1(i);
        u(2,i) = U2(i);
        xi(:,i)=dynamics(x,u(:,i),grid_world);
        J1s(i) = cost_to_go(xi(:,i),J);
        Js(i) = gamma*J1s(i) + cost_function(x);
    end
    [val, ind]= min(Js);
    u =u(:,ind);
end

%==========================================================================
% fit linear function (fill in)
%==========================================================================
function thetas = fit_linear(Js, xs, ys)
    states = [xs(:), ys(:)]';
    param = [];
    thetas = regularized_least_squares(@linear, states, Js, param);
end
 
%==========================================================================
% linear basis function (fill in)
%==========================================================================
function phi=linear(x,param)
    phi = [x(1); x(2); 1]; % [x, y, 1]
end

%==========================================================================
% quadratic basis function (fill in)
%==========================================================================
function phi=quadratic(x,param)
    phi = [x(1)^2; x(2)^2; x(1); x(2); x(1)*x(2); 1]; 
end

%==========================================================================
% fit quadratic function (fill in)
%==========================================================================
function thetas = fit_quadratic(Js, xs, ys)
    states = [xs(:), ys(:)]';
    param = [];
    thetas = regularized_least_squares(@quadratic, states, Js, param);
end

%==========================================================================
% regularized least squares
%==========================================================================
function theta=regularized_least_squares(BASIS,x,y,param)
    N=length(x);
    Phi=[];
    for k=1:N
        phi=BASIS(x(:,k),param);
        Phi=[Phi phi];
    end
    gamma=0.00;
    Phi=Phi';
    P=(Phi'*Phi);
    I=eye(length(P));
    B=Phi'*y;
    theta=(P+gamma*I)\B;
end

%==========================================================================
% fit cost-to-go
%==========================================================================
function Jc = fit_cost_to_go(Xi,Yi,Ji,Jc)
    switch Jc.type
        case 'linear'
            thetas = fit_linear(Ji, Xi, Yi);
        case 'quadratic'
            thetas = fit_quadratic(Ji, Xi, Yi);
    end    
    Jc.thetas = thetas;
end