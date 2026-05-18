%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements direct transcription with SNOPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,F,info] = dirtran(N,DYNAMICS,GRADIENTS,param)
global SNOPT_USERFUN;

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

A = [];
iAfun = [];
jAvar = [];

% Cost function bounds
Flow = -inf;
Fhigh = inf;
iGfun = ones(length(w0),1);
jGvar = [1:length(w0)]';
% f is the vector of costs + constraint evaluations f=[J;c1;c2...]
% G is the sparse matrix of cost and constraint gradients
% iGfun and iGvar represent the indices in G that are non-zero
nf = 1;

% set up the sparse structure
for i=1:N-1,
  for j=1:nX,
    nf = nf+1;
    iGfun = ; %FILL IN set up the spare structure for G
    jGvar = ;
    Flow = [Flow;0];
    Fhigh = [Fhigh;0];
  end
end

% set up the SNOPT parameters
setSNOPTParam('Major Iterations Limit',10000);
setSNOPTParam('Minor Iterations Limit',500);
setSNOPTParam('Major Optimality Tolerance',1e-6);
setSNOPTParam('Major Feasibility Tolerance',1e-6);
setSNOPTParam('Minor Feasibility Tolerance',1e-6);
setSNOPTParam('Superbasics Limit',200);
setSNOPTParam('Derivative Option',0);
setSNOPTParam('Verify Level',0);
setSNOPTParam('Iterations Limit',10000);

A = [];
iAfun = [];
jAvar = [];

%define the SNOPT userfunction
SNOPT_USERFUN = @(w) userFun(w,N,DYNAMICS,GRADIENTS,param);

%set up inputs for SNOPT solver
wmul = zeros(length(w0), 1);
wstate = zeros(length(w0), 1);
            
Fmul = zeros(length(Flow), 1);
Fstate = zeros(length(Flow), 1);

options.name = 'Prob0';
options.start = 'Cold';
options.screen = 'on';
options.printfile = ['prob0-', num2str(get(getCurrentTask(),'ID')), '.out'];

A_struct = struct('row', iAfun, 'col', jAvar, 'val', A);
G_struct = struct('row', iGfun, 'col', jGvar);
% call the SNOPT solver
[w, F, info] = snopt(w0, wlow, whigh, wmul, wstate, Flow, Fhigh, Fmul, Fstate, SNOPT_USERFUN, 0, 1, A_struct, G_struct, options);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for evaluating the cost, constraints, and gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,G] = userFun(w,N,DYNAMICS,GRADIENTS,param)
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
df_dyn = [];
f_u=[];
df_u=[];

J=0;
dJdx=[];
dJdu=[];

%construct the dynamics constraints and costs for each trajectory knot
%point
% f is the vector of costs + constraint evaluations f=[J;c1;c2...]
% G is the sparse matrix of cost and constraint gradients
for i=2:N
      
    t=0;
          
    xdot=DYNAMICS(t,X(:,i-1),u(:,i-1),param);
                   
    f_dyn(:,i-1) = ;%FILL IN the dynamics constraint (Assume EULER integration)
    
    [dfdx dfdu]=GRADIENTS(t,X(:,i-1),u(:,i-1),param);
    
    d_dyn=; % FILL IN the dynamics constraint gradients
           
    df_dyn(:,:,i-1) = d_dyn';
    
    [dJdxi dJdui]=; %FILL IN the running cost gradients
    
    Ji=cost(t,X(:,i-1),u(:,i-1),param); %evaluate the cost function
    
    J=J+Ji; %define the cost
    
    dJdx=[dJdx dJdxi]; %define the cost gradient w.r.t. state
    
    dJdu=[dJdu dJdui]; %define th cost gradient w.r.t. input

end
  

[dJdxi, dJdui]=final_cost_gradients(t,X(:,N),u(:,N),param); %compute the final cost gradients

Ji=final_cost(t,X(:,N),u(:,N),param); %compute the final cost

J=J+Ji;

dJdx=[dJdx dJdxi];% append the final cost gradients w.r.t. state

dJdu=[dJdu dJdui];%append the final cost gradients w.r.t input

cost0 = J; % the total cost
dcost = % FILL IN the cost gradients

f = [cost0;f_dyn(:)]; % cost and constraint evaluation
G = %FILL in G 
  
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the final cost gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dJdx dJdu]=final_cost_gradients(t,x,u,param)

Q=param.opt.Qf;
R=param.opt.R;

dJdx=2*x'*Q;
dJdu=0*2*u'*R;

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
% Compute the running cost gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dJdx dJdu]=cost_gradients(t,x,u,param)

Q=param.opt.Q;
R=param.opt.R;

dJdx=2*x'*Q;
dJdu=2*u'*R;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the running cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J=cost(t,x,u,param)

Q=param.opt.Q;
R=param.opt.R;

J=x'*Q*x+u'*R*u;

end


function setSNOPTParam(paramstring,default)
  snset([paramstring,'=',num2str(default)]);
end