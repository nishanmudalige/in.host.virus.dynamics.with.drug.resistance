omega = 0.05;
dV = 3;
rI = 0.0032;
lambda = 100;
dS = 0.1;
rR = 56.1;
R  = 1;
rP = 0.127;
P  = 1;
mR = 4.16;
mP = 8.52;
dI = 0.5;
nI = 62.5;
dR = 0; % 16.6
dP = 0; %8.32

% x(1) = VI
% x(2) = VNI
% x(3) = TS
% x(4) = TI
% x(5) = TR
% x(6) = TRP
% x(7) = TPNI
% x(8) = TPI


%for i = 1:10

    %epsilon = 0.1*i;

    g = @(t,x) [ ni*omega*x(4) - dv*x(1) - rI*x(3)*x(1) - rI*x(7)*x(1)  ; 
        nI*x(8) + nI*(1 - omega)*x(4) - dV*x(2) ;
        lambda - rI*x(3)*x(1) - dS*x(3) - rR*x(3)*R - rP*x(3)*P + mR*x(5) + mP*x(7) ;
        rI*x(3)*x(1) - dI*x(4) - rP*x(4)*P + mP*x(8) ;
        rR*x(3)*R - dS*x(5) + mP*x(6) - mR*x(5) - rP*x(6)*P ;
        rR*x(7)*R - dS*x(6) - mP*x(6) - mR*x(6) + rP*x(5)*P ;
        rP*x(3)*P - dS*x(7) - rI*x(7)*x(1) - rR*x(7)*R - mP*x(7) + mR*x(6) ;
        rI*x(7)*x(1) - dI*x(8) + rP*x(4)*P - mP*x(8) ;
        -dR*x(9);
        -dP*x(10) ] ;

    [t,xa] = ode45(g,[0 100],[100 0 1000 0 0 0 0 0 0 0]);

    % hold on
    plot(t, xa(:,1), t, xa(:,2), t, xa(:,3), t, xa(:,4), t, xa(:,5), t, xa(:,6), t, xa(:,7), t, xa(:,8) )
    
    
    %title(caption)
%legend('     Target', '     Sensitive', '     Resistant')
%legend('Location', 'NorthEast') % move legend to upper left

%end



%     g = @(t,x) [ ni*omega*x(4) - dv*x(1) - rI*x(3)*x(1) - rI*x(7)*x(1)  ; 
%         nI*x(8) + nI*(1 - omega)*x(4) - dV*x(2) ;
%         lambda - rI*x(3)*x(1) - dS*x(3) - rR*x(3)*R - rP*x(3)*P + mR*x(5) + mP*x(7) ;
%         rI*x(3)*x(1) - dI*x(4) - rP*x(4)*P + mP*x(8) ;
%         rR*x(3)*R - dS*x(5) + mP*x(6) - mR*x(5) - rP*x(6)*P ;
%         rR*x(7)*R - dS*x(6) - mP*x(6) - mR*x(6) + rP*x(5)*P ;
%         rP*x(3)*P - dS*x(7) - rI*x(7)*x(1) - rR*x(7)*R - mP*x(7) + mR*x(6) ;
%         rI*x(7)*x(1) - dI*x(8) + rP*x(4)*P - mP*x(8) ] ;
