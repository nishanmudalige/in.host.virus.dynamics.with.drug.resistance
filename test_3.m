% WIth drug and drug resistant virus

tn = 40*24; % hours
t0 = 0;  % initial time
h  = 0.001; 

t = t0:h:tn;

nt = length(t);

T   = zeros(nt,1);
TES = zeros(nt,1);
TS  = zeros(nt,1);
VS  = zeros(nt,1);
TER = zeros(nt,1);
TR  = zeros(nt,1);
VR  = zeros(nt,1);

% Initial conditions
T(1)   = 10^(6);
TES(1) = 0;
TS(1)  = 0;
VS(1)  = 10^(-6);
TER(1) = 0;
TR(1)  = 0;
VR(1)  = 0;

% Parameters
phiR = 0.8;
phiS = 1.25;
deltaE = 0.7;
deltaS = 1;
deltaRPR = deltaS;
b = 0.01;
k = 2.4*10^(-8);
lambda = 10^(4);
d = 0.01;
c = 23;
deltaE = 1;
ps = 4000;
pr = 4000;

sigmaRTI = 0.2;
sigmaPI = 0.2;

epsilonRTI = 0.3; %0.5;
epsilonPI = 0.3; %0.6;


% x(1) = target cells (T)
% x(2) = infected sensitive (Is)


for i = 1:nt-1
 
    T(i+1)   = T(i)   + h*( lambda - d*T(i) - k*VS(i)*T(i) - k*VR(i)*T(i) + b*TES(i) + b*TER(i) );
    TES(i+1) = TES(i) + h*( k*VS(i)*T(i) - (b + phiS*(1 - epsilonRTI) + deltaE)*TES(i) );
    TS(i+1)  = TS(i)  + h*( phiS*(1 - epsilonRTI)*TES(i) - deltaS*TS(i) );
    VS(i+1)  = VS(i)  + h*( ps*(1 - epsilonPI)*TS(i) - c*VS(i) );
    TER(i+1) = TER(i) + h*( k*VR(i)*T(i) - (b + phiR*(1 - epsilonRTI*sigmaRTI) + deltaE)*TER(i)  );
    TR(i+1)  = TR(i)  + h*( phiR*(1 - epsilonRTI*sigmaRTI)*TES(i) - deltaRPR*TR(i) );
    VR(i+1)  = VR(i)  + h*( pr*(1 - epsilonPI*sigmaPI)*TR(i) - c*VR(i) );

end

    
%plot(t,T, t, Is, t, Ir, 'LineWidth', 1.2)
%hold on

% plot(t, TS, t, TR, t, TES, 'Linewidth', 2)
% legend('  Resistant',...
%        '  Sensitive')
% legend('Location', 'NorthEast')

plot(t, VR, t, VS, 'Linewidth', 2)
legend('  Resistant',...
       '  Sensitive')
legend('Location', 'NorthEast')

