lambda = 5;
delta = 0.5;
b = 0.25;
a = 1;


% x(1) = target cells (T)
% x(2) = infected cells (I)

% format longG

for i = 1:10

    epsilon = 0.1*i;
    
    g = @(t,x) [lambda - delta*x(1) - (1-epsilon)*b*x(2)*x(1); (1-epsilon)*b*x(2)*x(1) - a*x(2) ];
    [t,xa] = ode45(g,[0 20],[lambda/delta 1]);
    
    subplot(2, 5, i)
    caption = sprintf('Epsilon = %d', epsilon); % PROBLEM !!!
    % caption = sprintf('Epsilon = %d', num2str(epsilon)); % PROBLEM !!!
    plot(t,xa(:,1), t, xa(:,2))
    title(caption)
    % hold on
    
    
end