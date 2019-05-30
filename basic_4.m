lambda = 5;
delta = 0.5;
bs = 0.25;
br = 0.24;
mu = 10^(-5);
a = 1;

% x(1) = target cells (T)
% x(2) = infected sensitive (Is)
% x(3) = infected resistant (Ir)


for i = 1:10
    
    %for j = 1:10

    epsilon1 = 0.1*i;
    epsilon2 = 0.8;%*j;

    g = @(t,x) [lambda - delta*x(1) - ((1-epsilon1)*(1-epsilon2)*bs*x(2) + br*x(3))*x(1) ; 
        (1-epsilon1)*(1-epsilon2)*(1-mu)*bs*x(2)*x(1) - a*x(2) ;
        br*x(1)*x(3) - a*x(3) + (1-epsilon1)*(1-epsilon2)*mu*bs*x(2)*x(1)] ;

    %[t,xa] = ode45(g,[0 20],[lambda/delta 1]);
    [t,xa] = ode45(g,[0 100],[4 3 0]);

    subplot(2, 5, i);
    caption = sprintf('E1 = %2g', epsilon1);
    plot(t,xa(:,1), 'k', t, xa(:,2), 'g', t, xa(:,3), 'r', 'LineWidth', 1.2)
    title(caption)
%legend('     Target', '     Sensitive', '     Resistant')
%legend('Location', 'NorthEast') % move legend to upper left

    %end

end