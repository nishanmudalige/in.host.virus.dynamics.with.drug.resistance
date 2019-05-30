% WIth timed drug
% clc
% clear

tn = 12*25; % hours
t0 = 0;  % initial time
h  = 0.001; 

t = t0:h:tn;

nt = length(t);

T    = zeros(nt,1);
IS   = zeros(nt,1);
IR1  = zeros(nt,1);
IR2  = zeros(nt,1);
IR12 = zeros(nt,1);

% Initial conditions
T(1)  = 4;
IS(1) = 3;
IR1(1) = 0;
IR2(1) = 0;
IR12(1) = 0;

% Parameters
lambda = 5/24;
delta = 0.5/24;
b = 0.25;
bs = 0.25;
br1 = 0.24;
br2 = 0.24;
br12 = 0.24;
mu = 10^(-5);
mu1 = (1/2)*mu;
mu2 = (1/3)*mu;
mu12 = (1/6)*mu;
a = 1;
%epsilon=0.9;

tvec = [2 4 8 12 24 48];

% time = 2;

% x(1) = target cells (T)
% x(2) = infected sensitive (Is)

    epsilon1 = 0.8;%*i;
    epsilon2 = 0.8;%*j;

    conc1 = 1.0;
    conc2 = 1.0;
    
%j = time/h;
%j = j/h;
k=0;
l=0;
m=0;

for j = 1:length(tvec)
    
    time = tvec(j)
    
    for i = 1:nt-1
    
        if( i < (time/h)  )

        k = k+1;

        conc1 = 0.1;
        conc2 = 0.1;

        T(i+1)  = T(i) + ...
            h*(lambda - delta*T(i) - ( conc1*(1-epsilon1)*conc2*(1-epsilon2)*bs*IS(i) + br1*IR1(i) + br2*IR2(i) + br12*IR12(i))*T(i) );
        IS(i+1) = IS(i) + ...
            h*( conc1*(1-epsilon1)*conc1*(1-epsilon2)*(1-mu)*bs*IS(i)*T(i) - a*IS(i) );
        IR1(i+1) = IR1(i) + ...
            h*( br1*T(i)*IR1(i) - a*IR1(i) + conc1*(1-epsilon1)*mu1*br1*IS(i)*T(i) );
        IR2(i+1) = IR2(i) + ...
            h*( br1*T(i)*IR2(i) - a*IR2(i) + conc2*(1-epsilon2)*mu2*br2*IS(i)*T(i) );
        IR12(i+1) = IR12(i) + ...
            h*( br1*T(i)*IR12(i) - a*IR12(i) + conc1*(1-epsilon1)*conc2*(1-epsilon2)*mu12*br12*IS(i)*T(i) );

        conc1 = 1.0;
        conc2 = 1.0;

        %elseif( mod((i/h),time) == 0 )
        %elseif( time * round(double( (i/h) )/time) == (i/h) )

        %elseif( mod(i,j) == 0) % & fix((i/h)/time) > time )
        %if( i == j )    
        elseif( mod(i,time/h) == 0 )

        l = l+1;

        conc1 = 2.5*10^(2);
        conc2 = 2.5*10^(2);

        T(i+1)  = T(i) + ...
            h*(lambda - delta*T(i) - ( conc1*(1-epsilon1)*conc2*(1-epsilon2)*bs*IS(i) + br1*IR1(i) + br2*IR2(i) + br12*IR12(i))*T(i) );
        IS(i+1) = IS(i) + ...
            h*( conc1*(1-epsilon1)*conc1*(1-epsilon2)*(1-mu)*bs*IS(i)*T(i) - a*IS(i) );
        IR1(i+1) = IR1(i) + ...
            h*( br1*T(i)*IR1(i) - a*IR1(i) + conc1*(1-epsilon1)*mu1*br1*IS(i)*T(i) );
        IR2(i+1) = IR2(i) + ...
            h*( br1*T(i)*IR2(i) - a*IR2(i) + conc2*(1-epsilon2)*mu2*br2*IS(i)*T(i) );
        IR12(i+1) = IR12(i) + ...
            h*( br1*T(i)*IR12(i) - a*IR12(i) + conc1*(1-epsilon1)*conc2*(1-epsilon2)*mu12*br12*IS(i)*T(i) );

        conc1 = 1.0;
        conc2 = 1.0;

        else

        m=m+1;

        conc1 = 1.0;
        conc2 = 1.0;

        T(i+1)  = T(i) + ...
            h*(lambda - delta*T(i) - ( conc1*(1-epsilon1)*conc2*(1-epsilon2)*bs*IS(i) + br1*IR1(i) + br2*IR2(i) + br12*IR12(i))*T(i) );
        IS(i+1) = IS(i) + ...
            h*( conc1*(1-epsilon1)*conc1*(1-epsilon2)*(1-mu)*bs*IS(i)*T(i) - a*IS(i) );
        IR1(i+1) = IR1(i) + ...
            h*( br1*T(i)*IR1(i) - a*IR1(i) + conc1*(1-epsilon1)*mu1*br1*IS(i)*T(i) );
        IR2(i+1) = IR2(i) + ...
            h*( br1*T(i)*IR2(i) - a*IR2(i) + conc2*(1-epsilon2)*mu2*br2*IS(i)*T(i) );
        IR12(i+1) = IR12(i) + ...
            h*( br1*T(i)*IR12(i) - a*IR12(i) + conc1*(1-epsilon1)*conc2*(1-epsilon2)*mu12*br12*IS(i)*T(i) );

        end % end if

    end % for i

    subplot(2, 3, j);

    plot(t,IR1, 'r',...
         t,IR2, 'c',...
         t,IR12, 'm',...
         'LineWidth', 1.2)
     xlabel('Time')
     ylabel('Infected')
     
     
   legend('  Resistant affected by drug 1',...
          '  Resistant affected by drug 2',...
          '  Resistant affected by drug 1 and 2')
   legend('Location', 'NorthEast')


end % for j

    %subplot(10, 10, i);
    %caption1 = sprintf('E1 = %2g', epsilon1);
    %caption2 = sprintf(' E2 = %2g', epsilon2);
    %plot(t,xa(:,1), 'k', t, xa(:,2), 'g', t, xa(:,3), 'r', 'LineWidth', 1.2)
    %title([caption1 caption2])
    
%     plot(t,T, 'b',...
%         t,IS, 'g',...

%    plot(t,IR1, 'r',...
%         t,IR2, 'c',...
%         t,IR12, 'm',...
%         'LineWidth', 1.2)
%    legend('  Resistant affected by drug 1',...
%           '  Resistant affected by drug 2',...
%           '  Resistant affected by drug 1 and 2')
%    legend('Location', 'NorthEast')
%    %
   %tn/time
   [fprintf('k ') fprintf('l ') fprintf('m ')]
   [k l m]
   
   fprintf('k+l+m=',k+l+m)
   (k+l+m) - tn/h
   %hold on

%    plot(t,T, ...
%           t,IS, 'LineWidth', 2)
%    legend('  Target',...
%           '  Sensitive')
%    legend('Location', 'NorthEast')

      
%     plot(t, IR1)
%     hold on
%     plot(t, IR2)
%     hold on
%     plot(t, IR12)
