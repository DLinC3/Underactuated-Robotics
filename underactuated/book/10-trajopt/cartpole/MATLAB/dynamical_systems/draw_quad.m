 function drawQuad(t,x,param)
      % Draw the quadrotor.  
      persistent hFig base pin prop;
  
      L=param.L;

      if (isempty(hFig))
        hFig = figure(25);
        set(hFig,'DoubleBuffer', 'on');
        
        base = [1.2*L*[1 -1 -1 1]; .025*[1 1 -1 -1]];
        pin = [.005*[1 1 -1 -1]; .1*[1 0 0 1]];
        a = linspace(0,2*pi,50);
        prop = [L/1.5*cos(a);.1+.02*sin(2*a)];
      end
            
      figure(hFig); cla; hold on; view(0,90);
      
      r = [cos(x(3)), -sin(x(3)); sin(x(3)), cos(x(3))];
      
      p = r*base;
      patch(x(1)+p(1,:), x(2)+p(2,:),1+0*p(1,:),'b','FaceColor',[.6 .6 .6])
      
      p = r*[L+pin(1,:);pin(2,:)];
      patch(x(1)+p(1,:),x(2)+p(2,:),0*p(1,:),'b','FaceColor',[0 0 0]);
      p = r*[-L+pin(1,:);pin(2,:)];
      patch(x(1)+p(1,:),x(2)+p(2,:),0*p(1,:),'b','FaceColor',[0 0 0]);
      
      p = r*[L+prop(1,:);prop(2,:)];
      patch(x(1)+p(1,:),x(2)+p(2,:),0*p(1,:),'b','FaceColor',[0 0 1]);
      p = r*[-L+prop(1,:);prop(2,:)];
      patch(x(1)+p(1,:),x(2)+p(2,:),0*p(1,:),'b','FaceColor',[0 0 1]);
      
      title(['t = ', num2str(t(1),'%.2f') ' sec']);
      %set(gca,'XTick',[],'YTick',[])
                     
      axis equal; axis([-1 1 -1 1]);
      title(['time: ',num2str(t,2),' s']);
      xlabel('x (m)'); ylabel('z (m)');
      
      
       
 end