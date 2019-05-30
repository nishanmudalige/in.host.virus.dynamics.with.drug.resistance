g = @(t,x) [lambda - delta*x(1) - (1-epsilon)*b*x(2)*x(1); (1-epsilon)*b*x(2)*x(1) - a*x(2) ]

syms T I lambda delta epsilon a b
jacobian([lambda - delta*T - (1-epsilon)*b*I*T; (1-epsilon)*b*I*T - a*I], [T, I])

% g = @(t,x) [lambda - delta*x(1) - (1-epsilon)*b*x(2)*x(1); (1-epsilon)*b*x(2)*x(1) - a*x(2) ]

% fx=diff(g,T)
% fy=diff(g,I)

% [xcr,ycr]=solve(fx,fy); [xcr,ycr]

fp = jacobian([lambda - delta*T - (1-epsilon)*b*I*T; (1-epsilon)*b*I*T - a*I], [T, I])
subs(fp, [T I], [a/(b*(1-epsilon) (k-delta)/(b*(1-epsilon)) ]