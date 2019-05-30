lambda = 5/24;
delta = 0.5/24;
bs = 0.25/24;
br = 0.24/24;
mu = 10^(-5)/24;
a = 1/24;

% x(1) = target cells (T)
% x(2) = infected sensitive (Is)
% x(3) = infected resistant (Ir)


for i = 1:(100)
    
    %for j = 1:5

    epsilon1 = 0.1;%*i;
    epsilon2 = 0.2;%*j;
    
    
    
    if( rem(i,6) == 0 )
    
    conc1 = 1.5;
    conc2 = 1.5;

    g = @(t,x) [lambda - delta*x(1) - (conc1*(1-epsilon1)*conc2*(1-epsilon2)*bs*x(2) + br*x(3))*x(1) ; 
        conc1*(1-epsilon2)*conc2*(1-epsilon2)*(1-mu)*bs*x(2)*x(1) - a*x(2) ;
        br*x(1)*x(3) - a*x(3) + conc1*(1-epsilon1)*conc2*(1-epsilon2)*mu*bs*x(2)*x(1)] ;

    else
        
    conc1 = 1.0;
    conc2 = 1.0;

    g = @(t,x) [lambda - delta*x(1) - (conc1*(1-epsilon1)*conc2*(1-epsilon2)*bs*x(2) + br*x(3))*x(1) ; 
        conc1*(1-epsilon2)*conc2*(1-epsilon2)*(1-mu)*bs*x(2)*x(1) - a*x(2) ;
        br*x(1)*x(3) - a*x(3) + conc1*(1-epsilon1)*conc2*(1-epsilon2)*mu*bs*x(2)*x(1)] ;
    end
    
    %[t,xa] = ode45(g,[0 20],[lambda/delta 1]);
    [t,xa] = ode45(g,[0 1000],[4 3 0]);

    subplot(10, 10, i);
    caption1 = sprintf('E1 = %2g', epsilon1);
    caption2 = sprintf(' E2 = %2g', epsilon2);
    plot(t,xa(:,1), 'k', t, xa(:,2), 'g', t, xa(:,3), 'r', 'LineWidth', 1.2)
    title([caption1 caption2])
%legend('     Target', '     Sensitive', '     Resistant')
%legend('Location', 'NorthEast') % move legend to upper left

    %end

end