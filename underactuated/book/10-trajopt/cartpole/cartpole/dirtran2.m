%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements direct transcription with FMINCON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,F,info] = dirtran2(N,DYNAMICS,param)

param.opt.xi=param.opt.x0;

nU=param.opt.nU;
nX=param.opt.nX;

%set up initial control input values
if isfield(param.opt,'u0s')
    u0=param.opt.u0s;
else
    u0 = zeros(nU*N,1);
end
% set up control input bounds
ulow = param.opt.umin*ones(nU*N,1);
uhigh = param.opt.umax*ones(nU*N,1);

x0=[];
xlow=[];
xhigh=[];

for k=1:nX   
    x0=[x0;linspace(param.opt.x0(k),param.opt.xd(k),N)]; %initial state values
    xlow=[xlow;param.opt.xmin(k)*ones(1,N)]; %state upper bounds
    xhigh=[xhigh;param.opt.xmax(k)*ones(1,N)]; %state lower bounds
end

if isfield(param.opt,'x0s')
    x0=param.opt.x0s;
else

end

%time step bounds
hhigh=param.opt.hmax;
hlow=param.opt.hmin;

xlow(:,1)=param.opt.x0;
xhigh(:,1)=param.opt.x0;
xlow(:,end)=param.opt.xfmin;
xhigh(:,end)=param.opt.xfmax;

% create single vector of optimization problem parameters
w0 = [2/N;x0(:);u0]; %initial parameter values
wlow = [hlow;xlow(:);ulow]; %parameter upper bounds
whigh = [hhigh;xhigh(:);uhigh]; %parameter lower bounds

% here add more tolerance to fmincon
options = optimset('Largescale','off','MaxFunEvals',300000,'Display','iter','TolFun',1e-3,'TolX',1e-3,'TolCon',1e-3);
info=1;
[w, F] = fmincon(@(w)cost_function(w,N,param), w0, [], [], [], [], wlow, whigh,@(w)constraints(w,N,DYNAMICS,param),options);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for evaluating the cost, constraints, and gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = cost_function(w,N,param)
persistent k;

if isempty(k)
    k=1;
end

nX=param.opt.nX;
nU=param.opt.nU;
I=eye(nX);

h = w(1);
X = reshape(w(2:(1+nX*N)),nX,N);
u = reshape(w(nX*N+2:end),nU,N);

f_dyn = [];

J=0;

for i=2:N
      
    t=0;          
    
    % h here denote dt
    J = J + cost(t, X(:,i-1), u(:,i-1), param) * h; %FILL IN the running cost (see helpful functions below)
    
end

% fprintf('running_cost = %f\n', J);
% since final_cost J=x'*Q*x;
% so x input should be x - xd

J = J + final_cost(t, X(:,N)-param.opt.xd, [], param); %FILL IN the total cost (see helpful functions below)

% fprintf('Terminal cost = %f\n', final_cost(t, X(:,N)-param.opt.xd, [], param));
cost0 = J; % the total cost


f = [cost0;f_dyn(:)]; % cost and constraint evaluation
  
if param.plot==1
    
    plot(0,0,'ko','MarkerSize',5,'MarkerFaceColor',[0 0 0]);
    xlabel('X (m)','FontSize',18);
    ylabel('Z (m)','FontSize',18);
    h=gca;
    set(h,'FontSize',16);
    axis equal; axis([-4 1 -1 1.75]);
    hold off    
    
end
 
k=k+1;
   
end

function [c,ceq] = constraints(w,N,DYNAMICS,param)

nX=param.opt.nX;
nU=param.opt.nU;
I=eye(nX);

h = w(1);
X = reshape(w(2:(1+nX*N)),nX,N);
u = reshape(w(nX*N+2:end),nU,N);

c =  [];
ceq= [];
for i=2:N
      
    t=0;
          
    xdot=DYNAMICS(t,X(:,i-1),u(:,i-1),param);

    f_dyn = X(:,i) - (X(:,i-1) + xdot * h); % FILL IN the dynamics constraint (Assume EULER integration)

    ceq = [ceq;f_dyn];
    
end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the final cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J=final_cost(t,x,u,param)

Q=param.opt.Qf;
R=param.opt.R;

J=x'*Q*x;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the running cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J=cost(t,x,u,param)

Q=param.opt.Q;
R=param.opt.R;

J=x'*Q*x+u'*R*u;

end