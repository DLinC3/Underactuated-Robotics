% ======================================================================
% This function implements the value iteration algorithm.
% It returns the optimal value function J and policy PI.
% ======================================================================

function [J,PI] = brick_vi_quadratic
global dt;

% define the mesh points as a row vector (YOU FILL THIS IN)
 q_bins = [-4:.1:4];   
 qdot_bins = [-5:.1:5];
if (isempty(q_bins) || isempty(qdot_bins)), error('you need to define the mesh'); end
tic;

% The discrete actions, defined as a row vector
a = 1*linspace(-1,1,9);
if (isempty(a)), error('you need to define the action set'); end

% dynamics dt
dt = 1e-2;

% create the mesh
[q qdot] = ndgrid(q_bins,qdot_bins);
s = [reshape(q,1,numel(q)); reshape(qdot,1,numel(qdot))];
ns = size(s,2); na = size(a,2);

% generate all possible state and action pairs
S = repmat(s,1,na); % repeat s na times
A = reshape(repmat(a,ns,1),1,ns*na); % repeat a ns times

% compute the one-step dynamics
Sn = S + dynamics(S,A).*dt;

% Compute the transition matrix
disp('Computing Transition Matrix...');
[Pi,P] = volumetric_interp(s,Sn,q_bins,qdot_bins);

% Compute the one-step cost
C = reshape(cost(S,A),ns,na);

% Setup value iteration (you shouldn't need to change this)
J = zeros(ns,1); % arbitrary initialization
converged  = 0.001; % value converged threshold
gamma = 0.9999; % discount factor
iter = 1; err = 1e6;

% Iterate the value estimate

% Value iteration update
while (err > converged)
    [Jnew, PI] = min(C + gamma*reshape(sum(P.*J(Pi),1),ns,na),[],2);
    err = max(abs(Jnew-J));
    disp(['iteration = ',num2str(iter),' ; max_err = ',num2str(err)]);
    J = Jnew; iter = iter+1;
    vi_plot_3d(J,a(PI),q_bins,qdot_bins,iter,'q','qdot');
end

disp('Value Estimate converged!');
% n1 = size(q_bins,2); n2 = size(qdot_bins,2);
% subplot(1,2,1);imagesc(q_bins,qdot_bins,reshape(PI,n1,n2)'); xlabel('q'); ylabel('q_{dot}'); title('u'); colorbar; axis xy; axis equal;
% subplot(1,2,2);imagesc(q_bins,qdot_bins,reshape(J,n1,n2)'); xlabel('q'); ylabel('q_{dot}'); title('Value function'); colorbar; axis xy; axis equal;
% 
% set(gcf, 'PaperPosition', [0 0 8 5]);%Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [8 5]); %Set the paper to have width 5 and height 5.
% % saveas(gcf, 'mintimevi', 'png') %Save figure

toc;
disp('Press Enter to simulate...');pause;
PI = a(PI)'; T = 10; disp_dts = 5;
xtraj = zeros(2,T/dt); x = [-1 1.2]';
for i=1:T/dt
    xtraj(:,i) = x;    
    t(i)=(i-1)*dt;
    [ind,coef] = volumetric_interp(s,x,q_bins,qdot_bins);
    u = sum(coef.*PI(ind),1);    
    utraj(i)=u;
    if (mod(i,disp_dts)==0)
        draw((i-1)*dt,x);
    end
    x = x + dynamics(x,u).*dt;
end



A = [0 1; 0 0];
B = [0; 1];

qp = 1;
qd = 1;
r = 10;

Q = diag([qp, qd]);  
R = r;               

[K_lqr, P_lqr] = lqr(A, B, Q, R);

% Calculate analytical cost-to-go J = x'*P*x
[q_mesh, qdot_mesh] = ndgrid(q_bins, qdot_bins);
J_lqr = zeros(size(q_mesh));
for i = 1:numel(q_mesh)
    x = [q_mesh(i); qdot_mesh(i)];
    J_lqr(i) = x'*P_lqr*x;
end

% Calculate analytical policy
PI_lqr = zeros(size(q_mesh));
for i = 1:numel(q_mesh)
    x = [q_mesh(i); qdot_mesh(i)];
    PI_lqr(i) = -K_lqr*x;
end

figure(1)
plot(t,xtraj(1,:),t,xtraj(2,:))
axis equal
xlabel('Time (s)','FontSize',14);
ylabel('System States','FontSize',14);
title('Minimum Time Policy Trajectory','FontSize',14);
legend(['Position','Velocity'])
hold off
set(gcf, 'PaperPosition', [0 0 8 5]);%Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, 'mintimetrajs', 'png') %Save figure


figure(2)
plot(xtraj(1,:),xtraj(2,:))
axis equal
xlabel('Position','FontSize',14);
ylabel('Velocity','FontSize',14);
title('Minimum Time Policy Trajectory','FontSize',14);
hold off
set(gcf, 'PaperPosition', [0 0 8 5]);%Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, 'mintimetraj', 'png') %Save figure

figure(3)
plot(t,utraj)
axis equal
xlabel('Time (s)','FontSize',14);
ylabel('Control Action','FontSize',14);
title('Minimum Time Control Action','FontSize',14);
hold off
set(gcf, 'PaperPosition', [0 0 8 5]);%Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, 'mintimecontrol', 'png') %Save figure

figure(4)
hold on
plot(-1,1.2,'r*');
plot(xtraj(1,:),xtraj(2,:),'--k')
axis equal
n1 = size(q_bins,2); n2 = size(qdot_bins,2);
imagesc(q_bins,qdot_bins,reshape(PI,n1,n2)'); axis xy;
ssurf_plot(q_bins);
hold off
xlabel('Position','FontSize',14);
ylabel('Velocity','FontSize',14);
title('Value Iteration vs. Analytical Solution','FontSize',14);
hold off
set(gcf, 'PaperPosition', [0 0 8 5]);%Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, 'mintimeanalytic', 'png') %Save figure
 
% ====================================================
% Plotting comparison
% ====================================================
figure(5)
subplot(2,2,1)
imagesc(q_bins,qdot_bins,reshape(PI,n1,n2)')
title('VI Policy'); axis xy; colorbar

subplot(2,2,2)
imagesc(q_bins,qdot_bins,reshape(J,n1,n2)')
title('VI Value Function'); axis xy; colorbar

% LQR Results
subplot(2,2,3)
imagesc(q_bins,qdot_bins,PI_lqr')
title('LQR Policy'); axis xy; colorbar

subplot(2,2,4)
imagesc(q_bins,qdot_bins,J_lqr')
title('LQR Value Function'); axis xy; colorbar

% Value function cross-section
figure(6)
plot(q_bins, J(1:length(q_bins)), 'b', q_bins, J_lqr(:,1), 'r--')
title('Value Function Comparison (qdot=0)')
legend('VI','LQR')
xlabel('Position'); ylabel('Value')


end % end of brick_vi

function ssurf_plot(q)

qdot=-sign(q).*sqrt(2*sign(q).*q);

plot(q,qdot,'k','LineWidth',4)


end

% ==========================================================
% This function performs volumetric interpolation on the
% state(s) Sn, returning the box indices(Pi) and weights(P)
% (a.k.a. transition probabilities)
% The input is:
% s: the mesh
% Sn: the states to interpolate for
% [q-bins, qdot_bins]: the bins used to create the mesh
% ==========================================================
function [Pi,P] = volumetric_interp(s,Sn,q_bins,qdot_bins)
ns = size(Sn,2);
Pi = zeros(4,ns);
P = Pi;

% impose limits and wrapping on the state
Sn = normalize(Sn,q_bins,qdot_bins);

% compute each transition individualy
for i=1:ns
    if (mod(i,4e3)==0)
        disp([num2str((i/ns)*1e2),'% done']);
    end
    % lower left corner
    ind_q = max([find(q_bins <= Sn(1,i),1,'last') 1]);
    ind_qdot = max([find(qdot_bins <= Sn(2,i),1,'last') 1]);
    
    offset = [0 0;1 0;0 1;1 1];
    if (ind_q == length(q_bins))
        offset(:,1) = -offset(:,1);
    end
    if (ind_qdot == length(qdot_bins))
        offset(:,2) = -offset(:,2);
    end
    % compute the total area of the containing box
    totl_area = abs(q_bins(ind_q+offset(2,1))-q_bins(ind_q))*...
                abs(qdot_bins(ind_qdot+offset(3,2))-qdot_bins(ind_qdot));
    % compute the four corner indices and weights
    for j=1:4
        state = [q_bins(ind_q+offset(j,1)); qdot_bins(ind_qdot+offset(j,2))];
        Pi(j,i) = find(sum(abs(s-repmat(state,1,size(s,2))),1)==0);
        P(5-j,i) = (abs(Sn(1,i)-q_bins(ind_q+offset(j,1)))*abs(Sn(2,i)-qdot_bins(ind_qdot+offset(j,2))))/totl_area;
    end
end
end

% =============================================================
% This function defines the continuous dynamics
% If you change dynamics, make sure that it is in vector form
% CHANGE DYNAMICS FOR PART (c)
% =============================================================
function xdot = dynamics(x,u)
xdot = [x(2,:); u];
end

% ===================================================================h
% This function defines the instantaneous cost (i.e. g(x,u))
% YOU SHOULD FILL THIS IN. 
% Note that X and u are vectors
% ===================================================================
function C = cost(X,u)
global dt;
qp = 1; 
qd = 1; 
r = 10;
C = 0.5 * ( qp*(X(1,:).^2) + qd*(X(2,:).^2) + r*(u.^2) ) * dt;
end

% ==============================================================
% This function imposes limits on the state
% ==============================================================
function s = normalize(s,q_bins,qdot_bins)
% impose limits
N = size(s,2);
smax = repmat([q_bins(end);qdot_bins(end)],1,N); ind = s>smax;
s(ind) = smax(ind);
smin = repmat([q_bins(1);qdot_bins(1)],1,N); ind = s<smin;
s(ind) = smin(ind);
end

% ===============================================================
% This function plots the value function and policy in 3D
% during value iteration.
% ===============================================================
function vi_plot_3d(J,PI,q_bins,qdot_bins,iter,xlabel_name,ylabel_name)

figure(10);
clf;

n1 = length(q_bins);
n2 = length(qdot_bins);

% The original code stores J and PI using ndgrid(q_bins,qdot_bins),
% so reshape(...,n1,n2)' gives matrices whose rows correspond to qdot
% and columns correspond to q.
J_surf  = reshape(J,n1,n2)';
PI_surf = reshape(PI,n1,n2)';

[Q,Qdot] = meshgrid(q_bins,qdot_bins);

% -----------------------------
% Cost-to-go surface
% -----------------------------
subplot(1,2,1);
surf(Q,Qdot,J_surf);
shading interp;
xlabel(xlabel_name);
ylabel(ylabel_name);
zlabel('J');
title('Cost-to-Go');
axis tight;
grid on;
view(-35,30);

% -----------------------------
% Policy surface
% -----------------------------
subplot(1,2,2);
surf(Q,Qdot,PI_surf);
shading interp;
xlabel(xlabel_name);
ylabel(ylabel_name);
zlabel('u');
title(['Policy, iteration ',num2str(iter)]);
axis tight;
grid on;
view(-35,30);

drawnow;

end

% ==============================================================
% This is the draw function for the brick
% ==============================================================
function draw(t,x)

persistent hFig blockx blocky;

if (isempty(hFig))
  hFig = figure(25);
  set(hFig,'DoubleBuffer','on');
  blockx = [-1, -1, 1, 1, -1];
  blocky = [0, 0.5, 0.5, 0, 0];
end

figure(hFig);
clf;

% draw the mass
brickcolor=[.75 .6 .5];
fill(blockx+repmat(x(1),1,5),blocky,brickcolor);
hold on

faintline=[.6 .8 .65]*1.1;
plot(min(blockx)+[0 0],[-5 5],'k:','Color',faintline);
plot(max(blockx)+[0 0],[-5 5],'k:','Color',faintline);

% draw the ground
line([-5, 5], [0, 0],'Color',[.3 .5 1],'LineWidth',1);
axis([-5 5 -1 2]);
%grid on
axis equal;
title(['t = ', num2str(t)]);

drawnow;
end
