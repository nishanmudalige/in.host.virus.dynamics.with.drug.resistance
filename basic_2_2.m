lambda = 5;
delta = 0.5;
bs = 0.25;
br = 0.20;
mu = 10^(-5);
a = 1;

% x(1) = target cells (T)
% x(2) = infected sensitive (Is)
% x(3) = infected resistant (Ir)

%set(gca, 'XTickLabel', [],'XTick',[])



for i = 0:5

    epsilon = 0.2*i;

    g = @(t,x) [lambda - delta*x(1) - ((1-epsilon)*bs*x(2) + br*x(3))*x(1) ; 
        (1-epsilon)*(1-mu)*bs*x(2)*x(1) - a*x(2) ;
        br*x(1)*x(3) - a*x(3) + (1-epsilon)*mu*bs*x(2)*x(1)] ;

    %[t,xa] = ode45(g,[0 20],[lambda/delta 1]);
    [t,xa] = ode45(g,[0 100],[4 3 0]);

    subplot(2, 3, i+1);
    caption = sprintf('Epsilon = %2g', epsilon);    
    
    
    plot(t,xa(:,1), 'b', 'LineWidth', 1.2)
    hold on
    plot(t, xa(:,2), 'Color', [0 .5 0], 'LineWidth', 1.2)
    hold on
    plot(t, xa(:,3), 'r', 'LineWidth', 1.2)
    
    % hold on
    % total;
    
    
    title(caption)
    legend('     Target', '     Sensitive', '     Resistant')
    legend('Location', 'SouthOutSide') % move legend to upper left

end

