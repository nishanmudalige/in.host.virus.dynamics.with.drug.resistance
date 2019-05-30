lambda = 5;
delta = 0.5;
bs = 0.25;
br1 = 0.24;
br2 = 0.24;
br12 = 0.24;
mu = 10^(-5);
mu1 = (1/2)*mu;
mu2 = (1/3)*mu;
mu12 = (1/6)*mu;
a = 1;

% x(1) = target cells (T)
% x(2) = infected sensitive (IS)
% x(3) = infected resistant affected only by drug 1 (IR1)
% x(4) = infected resistant affected only by drug 2 (IR2)
% x(5) = infected resistant affected only by drug 1 and 2 (IR12)

%for i = 1:10
    
    %for j = 1:10

    epsilon1 = 1.0;%*i;
    epsilon2 = 1.0;%*j;
    
    g = @(t,x) [lambda - delta*x(1) - ((1-epsilon1)*(1-epsilon2)*bs*x(2) + br1*x(3) + br2*x(4) + br12*x(5))*x(1) ;
       (1-epsilon1)*(1-epsilon2)*(1-mu)*bs*x(2)*x(1) - a*x(2) ;
       br1*x(1)*x(3) - a*x(3) + (1-epsilon1)*mu1*br1*x(2)*x(1);
       br2*x(1)*x(4) - a*x(4) + (1-epsilon2)*mu2*br2*x(2)*x(1);
       br12*x(1)*x(5) - a*x(5) + (1-epsilon1)*(1-epsilon2)*mu12*br12*x(2)*x(1);
       ] ;

    %[t,xa] = ode45(g,[0 20],[lambda/delta 1]);
   [t,xa] = ode45(g,[0 25],[10 1 0 0 0]);

    %subplot(2, 5, i);
   % caption = sprintf('E1 = %2g', epsilon1);
   %plot(t,xa(:,1), 'k', t, xa(:,2), 'g', t, xa(:,3), 'r', 'LineWidth', 1.2)
   plot(t,xa(:,1), 'b', 'LineWidth', 2)
   hold on
   plot(t,xa(:,2), 'Color', [1 .5 0], 'LineWidth', 2)
   hold on
   plot(t,xa(:,3), 'g', 'LineWidth', 2)
   hold on
   plot(t,xa(:,4), 'm', 'LineWidth', 2)
   hold on
   plot(t,xa(:,5), 'c', 'LineWidth', 2)
   
   % title(caption)
legend('  Target', '  Sensitive',...
    '  Resistant affected by drug 1',...
    '  Resistant affected by drug 2',...
    '  Resistant affected by drug 1 and 2')
legend('Location', 'NorthEast') % move legend to upper left

    %end

%end




%     g = @(t,x) [lambda - delta*T - ((1-epsilon1)*(1-epsilon2)*bs*Is 
%                 + br1*Ir1 + br2*Ir2 + br12*Ir12)*T ; 
%         
%        (1-epsilon1)*(1-epsilon2)*(1-mu)*bs*IS*T - a*IS ;
%        
%        br1*T*IR1 - a*IR1 + (1-epsilon1)*mu1*bs*IS*T;
%        
%        br2*T*IR2 - a*IR2 + (1-epsilon2)*mu2*bs*IS*T;
%        
%        br12*T*IR12 - a*IR12 + (1-epsilon1)*(1-epsilon2)*mu12*bs*IS*T;
%        ] ;
