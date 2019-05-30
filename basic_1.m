%clc
%clear all
%format long

lambda = 5;
delta = 0.5;
b = 0.25;
a = 1;


% x(1) = target cells (T)
% x(2) = infected cells (I)

% format longG

PaperSize = [11 8.5];

for i = 0:5

    epsilon = 0.2*i;
    
    
    
    % T I
    g = @(t,x) [lambda - delta*x(1) - (1-epsilon)*b*x(2)*x(1); (1-epsilon)*b*x(2)*x(1) - a*x(2) ];
    [t,xa] = ode45(g,[0 30],[10 1]);
    
    %Disease free
    %[t,xa] = ode45(g,[0 30],[(lambda/delta) 0]);
    
    %Endemic
    %[t,xa] = ode45(g,[0 30],[-a/(b*(-1+epsilon)) (-lambda*b+lambda*b*epsilon+delta*a)/(a*b*(-1+epsilon))]);
    
    %figure;
    subplot(2, 3, i+1);
    caption = sprintf('Epsilon = %2g', epsilon);
    % caption = sprintf('Epsilon = %d', epsilon); % PROBLEM !!!
    % caption = sprintf('Epsilon = %d', num2str(epsilon)); % PROBLEM !!!
    plot(t,xa(:,1), t, xa(:,2), 'LineWidth', 2)
    title(caption)
    
    legend('     Target', '     Infected')
    legend('Location', 'NorthEast') % move legend to upper left
    
    %hold on
    
    
end

