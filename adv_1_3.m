% WIth timed drug

tn = 20*24; % hours
t0 = 0;  % initial time
h  = 0.001; 

t = t0:h:tn;

nt = length(t);

T = zeros(nt,1);
Is = zeros(nt,1);
Ir = zeros(nt,1);

% Initial conditions
T(1)  = 4;
Is(1) = 3;
Ir(1) = 0;

% Parameters
lambda = 5/24;
delta = 0.5/24;
b = 0.25/24;
bs = 0.25/24;
br = 0.24/24;
mu = 10^(-5)/24;
a = 1/24;
epsilon=0.9;

% x(1) = target cells (T)
% x(2) = infected sensitive (Is)


for i = 1:nt-1
    

    if( rem(i,4) == 0 )

    conc = 10.0;
        
    T(i+1)  = T(i) + h*(lambda - delta*T(i) - (conc*(1-epsilon)*bs*Is(i) + br*Ir(i))*T(i) );
    Is(i+1) = Is(i) + h*( conc*(1-epsilon)*(1-mu)*bs*Is(i)*T(i) - a*Is(i) );
    Ir(i+1) = Ir(i) + h*( br*T(i)*Ir(i) - a*Ir(i) + conc*(1-epsilon)*mu*bs*Is(i)*T(i) );
        
    else
    
    T(i+1)  = T(i) + h*(lambda - delta*T(i) - ((1-epsilon)*bs*Is(i) + br*Ir(i))*T(i) );
    Is(i+1) = Is(i) + h*( (1-epsilon)*(1-mu)*bs*Is(i)*T(i) - a*Is(i) );
    Ir(i+1) = Ir(i) + h*( br*T(i)*Ir(i) - a*Ir(i) + (1-epsilon)*mu*bs*Is(i)*T(i) );
        
    end

end

    %subplot(10, 10, i);
    %caption1 = sprintf('E1 = %2g', epsilon1);
    %caption2 = sprintf(' E2 = %2g', epsilon2);
    %plot(t,xa(:,1), 'k', t, xa(:,2), 'g', t, xa(:,3), 'r', 'LineWidth', 1.2)
    %title([caption1 caption2])
    
    plot(t,T, t, Is, t, Ir, 'LineWidth', 1.2)
    hold on
