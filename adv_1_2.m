% WIth timed drug

tn = 10*24; % hours
t0 = 0;  % initial time
h  = 0.001; 

t = t0:h:tn;

nt = length(t);

T = zeros(nt,1);
I = zeros(nt,1);

% Initial conditions
T(1) = lambda/delta;
I(1) = 1;

% Parameters
lambda = 5/24;
delta = 0.5/24;
b = 0.25/24;
bs = 0.25/24;
br = 0.24/24;
mu = 10^(-5)/24;
a = 1/24;
epsilon=0.1;

% x(1) = target cells (T)
% x(2) = infected sensitive (Is)


for i = 1:nt-1
    

    if( rem(i,6) == 0 )

    conc = 1.0;
        
    T(i+1) = T(i) + h*(lambda - delta*T(i) - conc*(1-epsilon)*b*I(i)*T(i) );
    I(i+1) = I(i) + h*( conc*(1-epsilon)*b*I(i)*T(i) - a*I(i) );

        
    else
    
    T(i+1) = T(i) + h*(lambda - delta*T(i) - (1-epsilon)*b*I(i)*T(i) );
    I(i+1) = I(i) + h*( (1-epsilon)*b*I(i)*T(i) - a*I(i) );

        
    end

end

    %subplot(10, 10, i);
    %caption1 = sprintf('E1 = %2g', epsilon1);
    %caption2 = sprintf(' E2 = %2g', epsilon2);
    %plot(t,xa(:,1), 'k', t, xa(:,2), 'g', t, xa(:,3), 'r', 'LineWidth', 1.2)
    %title([caption1 caption2])
    
    plot(t,T, t, I)
