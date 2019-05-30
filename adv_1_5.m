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

epsilon1=0.5;
epsilon2=0.5;

% x(1) = target cells (T)
% x(2) = infected sensitive (Is)

time = 4;

for i = 1:nt-1
    
    if (i>=0 & i<= time)
    
    conc1 = 0.0;
    conc2 = 0.0;
        
    T(i+1)  = T(i) + h*(lambda - delta*T(i) - (conc1*(1-epsilon1)*conc2*(1-epsilon2)*bs*Is(i) + br*Ir(i))*T(i) );
    Is(i+1) = Is(i) + h*( conc1*(1-epsilon1)*conc2*(1-epsilon2)*(1-mu)*bs*Is(i)*T(i) - a*Is(i) );
    Ir(i+1) = Ir(i) + h*( br*T(i)*Ir(i) - a*Ir(i) + conc1*(1-epsilon1)*conc2*(1-epsilon2)*mu*bs*Is(i)*T(i) );


    else if( rem(i,time) == 0 )

    conc1 = 2.0;
    conc2 = 2.0;
        
    T(i+1)  = T(i) + h*(lambda - delta*T(i) - (conc1*(1-epsilon1)*conc2*(1-epsilon2)*bs*Is(i) + br*Ir(i))*T(i) );
    Is(i+1) = Is(i) + h*( conc1*(1-epsilon1)*conc2*(1-epsilon2)*(1-mu)*bs*Is(i)*T(i) - a*Is(i) );
    Ir(i+1) = Ir(i) + h*( br*T(i)*Ir(i) - a*Ir(i) + conc1*(1-epsilon1)*conc2*(1-epsilon2)*mu*bs*Is(i)*T(i) );
        
    else
    
    T(i+1)  = T(i) + h*(lambda - delta*T(i) - ((1-epsilon1)*(1-epsilon2)*bs*Is(i) + br*Ir(i))*T(i) );
    Is(i+1) = Is(i) + h*( (1-epsilon1)*(1-epsilon2)*(1-mu)*bs*Is(i)*T(i) - a*Is(i) );
    Ir(i+1) = Ir(i) + h*( br*T(i)*Ir(i) - a*Ir(i) + (1-epsilon1)*(1-epsilon2)*mu*bs*Is(i)*T(i) );
        
    end

end

    %subplot(10, 10, i);
    %caption1 = sprintf('E1 = %2g', epsilon1);
    %caption2 = sprintf(' E2 = %2g', epsilon2);
    %plot(t,xa(:,1), 'k', t, xa(:,2), 'g', t, xa(:,3), 'r', 'LineWidth', 1.2)
    %title([caption1 caption2])
    
    plot(t,T, t, Is, t, Ir, 'LineWidth', 1.2)
    hold on
