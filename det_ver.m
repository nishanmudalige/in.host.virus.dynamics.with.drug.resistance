%clc
%clear all
%format long

betaS = 2.8*10^(-5) ;
betaV = 0.5*betaS ;
betaR = 0.5*betaS ;
d = 4 ;
c = 3 ;
thetaS = 1.2*10^(-2) ;
thetaR = 1.2*10^(-2) ;
p = 0.9 ;


% x(1) = T
% x(2) = IS
% x(3) = IR
% x(4) = VS
% x(5) = VR

% format longG

%for i = 1:10

%   epsilon = 0.1*i;
    
    g = @(t,x) [-betaS*x(4)*x(1) - betaV*x(5)*x(1) ;
         betaS*x(4)*x(1) - d*x(2) ;
         betaR*x(5)*x(1) - d*x(3) ; 
         p*thetaS*x(2) - c*x(4) ;
         (1-p)*thetaS*x(2) + thetaR*x(3) - c*x(5)];
     
    %[t,xa] = ode45(g,[0 30],[4*10^(8) 0 0  9.3*10^(-2) ]);
    [t,xa] = ode45(g,[0 30], [(4*10^8) 100 100  (9.3*10^(-2)) 100]);
    %[t,xa] = ode45(g,[0 30], [1 1 1 1 1]);
    
    %subplot(2, 5, i);
    %caption = sprintf('Epsilon = %2g', epsilon);
    % caption = sprintf('Epsilon = %d', epsilon); % PROBLEM !!!
    % caption = sprintf('Epsilon = %d', num2str(epsilon)); % PROBLEM !!!
    plot(t,xa(:,1), t, xa(:,2), t, xa(:,3), t, xa(:,4), t, xa(:,5) )
    %title(caption)
    % hold on
    
    
%end