function dydt = react(t,y)

% Solve the kinetics example
dydt = zeros(size(y));

% Parameters - reaction-rate constants
k1 = 5; k2 = 2; k3 = 1;
A = y(1); We'll be explicit about it here though you can do
B = y(2); the calculations directly with the y-values.
C = y(3);

% Evaluate the RHS expression
dydt(1) = -k1*A + k2*B;
dydt(2) = k1*A - (k2+k3)*B;
dydt(3) = k3*B;

