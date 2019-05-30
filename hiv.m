lambda = 5;
delta = 0.5;
b = 0.25;
a = 1;
epsilon = 0.1;

% x(1) = target cells (T)
% x(2) = infected cells (I)

g = @(t,x) [lambda - delta*x(1) - (1-epsilon)*b*x(2)*x(1); (1-epsilon)*b*x(2)*x(1) - a*x(2) ];
    [t,xa] = ode45(g,[0 20],[lambda/delta 1]);

    
plot(t,xa(:,1), t, xa(:,2))