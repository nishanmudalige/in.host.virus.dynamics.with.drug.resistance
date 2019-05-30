% WIth timed drug

tn = 20*24; % hours
t0 = 0;  % initial time
h  = 0.001; 

t = t0:h:tn;

nt = length(t);

T  = zeros(nt,1);
TE = zeros(nt,1);
I  = zeros(nt,1);
V  = zeros(nt,1);


% Initial conditions
T(1)  = 10^(6);
TE(1) = 0;
I(1) = 0;
V(1) = 10^(-6);

% Parameters
phi = 0.8;
deltaE = 0.7;
b = 0.01;
k = 2.4*10^(-8);
lambda = 10^(4);
d = 0.01;
c = 23;
delta = 1;
p = 4000;

% x(1) = target cells (T)
% x(2) = infected sensitive (Is)


for i = 1:nt-1
    
 
    T(i+1)  = T(i)  + h*( lambda - d*T(i) - k*V(i)*T(i) + b*TE(i) );
    TE(i+1) = TE(i) + h*( k*V(i)*T(i) - (b + phi + delta)*TE(i) );
    I(i+1)  = I(i)  + h*( phi*TE(i) - delta*I(i) );
    V(i+1)  = V(i)  + h*( p*I(i) - c*V(i) );

end

    
%plot(t,T, t, Is, t, Ir, 'LineWidth', 1.2)
%hold on

plot(t, V)
