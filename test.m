tn = 20; % hours
t0 = 0;  % initial time
h  = 0.001; 

t = t0:h:tn;

nt = length(t);

VI = zeros(nt,1);
VNI = zeros(nt,1);
TS = zeros(nt,1);
TI = zeros(nt,1);
TR = zeros(nt,1);
TRP = zeros(nt,1);
TPNI = zeros(nt,1);
TPI = zeros(nt,1);

% drugs
R = zeros(nt,1);
P = zeros(nt,1);

% Initial conditions
VI(1)   = 10;
VNI(1)  = 1;
TS(1)   = 10000;
TI(1)   = 1;
TR(1)   = 1;
TRP(1)  = 1;
TPNI(1) = 1;
TPI(1)  = 1;

R(1) = 0;
P(1) = 0;

omega = 0.05;
dV = 3;
rI = 0.0032;
lambda = 100;
dS = 0.1;
rR = 56.1;
R  = 1;
rP = 0.127;
P  = 1;
mR = 4.16;
mP = 8.52;
dI = 0.5;
nI = 62.5;

dR = 16.6;
dP = 8.32;
Ri = 7.3;
Pi = 11.6;

% VI(i) = VI
% VNI(i) = VNI
% TS(i) = TS
% TI(i) = TI
% TR(i) = TR
% TRP(i) = TRP
% TPNI(i) = TPNI
% TPI(i) = TPI


%for i = 1:10

    %epsilon = 0.1*i;

%     g = @(t,x) [ ni*omega*TI(i) - dv*VI(i) - rI*TS(i)*VI(i) - rI*TPNI(i)*VI(i)  ; 
%         nI*TPI(i) + nI*(1 - omega)*TI(i) - dV*VNI(i) ;
%         lambda - rI*TS(i)*VI(i) - dS*TS(i) - rR*TS(i)*R - rP*TS(i)*P + mR*TR(i) + mP*TPNI(i) ;
%         rI*TS(i)*VI(i) - dI*TI(i) - rP*TI(i)*P + mP*TPI(i) ;
%         rR*TS(i)*R - dS*TR(i) + mP*TRP(i) - mR*TR(i) - rP*TRP(i)*P ;
%         rR*TPNI(i)*R - dS*TRP(i) - mP*TRP(i) - mR*TRP(i) + rP*TR(i)*P ;
%         rP*TS(i)*P - dS*TPNI(i) - rI*TPNI(i)*VI(i) - rR*TPNI(i)*R - mP*TPNI(i) + mR*TRP(i) ;
%         rI*TPNI(i)*VI(i) - dI*TPI(i) + rP*TI(i)*P - mP*TPI(i) ;
%         -dR*x(9);
%         -dP*x(10) ] ;
% 
%     [t,xa] = ode45(g,[0 100],[100 0 1000 0 0 0 0 0 0 0]);

for i = 1:nt-1
    

    if( rem(i, 2/h) == 0 )

    VI(i+1) = VI(i) + h*( ni*omega*TI(i) - dv*VI(i) - rI*TS(i)*VI(i) - rI*TPNI(i)*VI(i) );
    VNI(i+1) = VNI(i) + h*( nI*TPI(i) + nI*(1 - omega)*TI(i) - dV*VNI(i) );
    TS(i+1) = TS(i) + h*( lambda - rI*TS(i)*VI(i) - dS*TS(i) - rR*TS(i)*R(i) - rP*TS(i)*P(i) + mR*TR(i) + mP*TPNI(i) );
    TI(i+1) = TI(i) + h*( rI*TS(i)*VI(i) - dI*TI(i) - rP*TI(i)*P(i) + mP*TPI(i) );
    TR(i+1) = TR(i) + h*( rR*TS(i)*R(i) - dS*TR(i) + mP*TRP(i) - mR*TR(i) - rP*TRP(i)*P(i) );
    TRP(i+1) = TRP(i) + h*( rR*TPNI(i)*R(i) - dS*TRP(i) - mP*TRP(i) - mR*TRP(i) + rP*TR(i)*P(i) );
    TPNI(i+1) = TPNI(i) + h*( rP*TS(i)*P(i) - dS*TPNI(i) - rI*TPNI(i)*VI(i) - rR*TPNI(i)*R(i) - mP*TPNI(i) + mR*TRP(i) );
    TPI(i+1) = TPI(i) + h*( rI*TPNI(i)*VI(i) - dI*TPI(i) + rP*TI(i)*(i) - mP*TPI(i) );
    
    R(i+1) = R(i) + Ri ;
    P(i+1) = P(i) + Pi ;
    
    else
    
    VI(i+1) = VI(i) + h*( ni*omega*TI(i) - dv*VI(i) - rI*TS(i)*VI(i) - rI*TPNI(i)*VI(i) );
    VNI(i+1) = VNI(i) + h*( nI*TPI(i) + nI*(1 - omega)*TI(i) - dV*VNI(i) );
    TS(i+1) = TS(i) + h*( lambda - rI*TS(i)*VI(i) - dS*TS(i) - rR*TS(i)*R(i) - rP*TS(i)*P(i) + mR*TR(i) + mP*TPNI(i) );
    TI(i+1) = TI(i) + h*( rI*TS(i)*VI(i) - dI*TI(i) - rP*TI(i)*P(i) + mP*TPI(i) );
    TR(i+1) = TR(i) + h*( rR*TS(i)*R(i) - dS*TR(i) + mP*TRP(i) - mR*TR(i) - rP*TRP(i)*P(i) );
    TRP(i+1) = TRP(i) + h*( rR*TPNI(i)*R(i) - dS*TRP(i) - mP*TRP(i) - mR*TRP(i) + rP*TR(i)*P(i) );
    TPNI(i+1) = TPNI(i) + h*( rP*TS(i)*P(i) - dS*TPNI(i) - rI*TPNI(i)*VI(i) - rR*TPNI(i)*R(i) - mP*TPNI(i) + mR*TRP(i) );
    TPI(i+1) = TPI(i) + h*( rI*TPNI(i)*VI(i) - dI*TPI(i) + rP*TI(i)*(i) - mP*TPI(i) );
    
    R(i+1) = R(i) + h*( -dR*R(i) ) ;
    P(i+1) = P(i) + h*( -dP*P(i) ) ;
    
     end

end
    
    % hold on
   
plot(t, TS)
    
    
    %title(caption)
%legend('     Target', '     Sensitive', '     Resistant')
%legend('Location', 'NorthEast') % move legend to upper left

%end



