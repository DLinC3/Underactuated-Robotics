% ============================================
% Draw function for the glider
% ============================================
function draw_glider(t,x,param)
xd=param.xd;

% x: [x z theta phi dot(x,z,theta,phi)]
sc = 2;
lw = -0.15*sc; %lw = -0.03*sc;
lh = 0.45*sc;
le = 0.04*sc;
mac = .1145*sc;

ct = cos(x(3)); st = sin(x(3));
ctp = cos(x(3)+x(4)); stp = sin(x(3)+x(4));

figure(25); clf; hold on;

%cg
plot(x(1),x(2),'r.','MarkerSize',10);

% fuselage
if (lw < 0)
    line([x(1) x(1)-lh*ct],[x(2) x(2)-lh*st],...
        'LineWidth',1,'Color',[0 0 0]);
else
    line([x(1)+(lw+mac/2)*ct x(1)-lh*ct],...
        [x(2)+(lw+mac/2)*st x(2)-lh*st],...
        'LineWidth',1,'Color',[0 0 0]);
end

% wing
line([x(1)+(lw+mac/2)*ct x(1)+(lw-mac/2)*ct],...
    [x(2)+(lw+mac/2)*st x(2)+(lw-mac/2)*st],...
    'LineWidth',3,'Color',[0 .2 1]);

% elevator
line([x(1)-lh*ct x(1)-lh*ct-2*le*ctp],...
    [x(2)-lh*st x(2)-lh*st-2*le*stp],...
    'LineWidth',3,'Color',[0 .2 1]);

axis equal; axis([-4 1 -1 1.75]);
title(['time: ',num2str(t,2),' s']);
xlabel('x (m)'); ylabel('z (m)');

% the perch
plot(xd(1),xd(2),'ko','MarkerSize',5,'MarkerFaceColor',[0 0 0]);
end