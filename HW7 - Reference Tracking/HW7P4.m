% ------------------------------------------------------------------------ 
%
%  ECE171A: HW7 Problem 4  -- sample code
%  Writen by Yang Zheng
%
% ------------------------------------------------------------------------ 

% Note: you can also simulate this problem using ODE
%       as what you have done in previous assignments.
%       For LTI systems, we can do the simulations using
%       transfer functions which are easier for coding. 

close all
%% 4b)
% System parameters
m = 1000;   % mass
c = 50;     % damping coefficient
b = 25;     % conversion factor 
a = 0.2;    % lag coefficient
T = 200;    % another conversion factor

% Dynamics 
s = tf('s');

P1 = b/(m*s+c);  % transfer function for tau to v
P2 = a*T/(s+a);  % transfer function for u to tau

%Case 1: kp=0.01:
% P controller 
kp  = 0.01; 
ki  = 0; %0.005;
C   = kp + ki/s;

L   = C*P2*P1;        % Loop transfer function 
Gyr1 = feedback(L,1);  % closed transfer function from r to y

% step response
figure; 
step(Gyr1);
title('Step Response (k_p=0.01)');

% freqency response 
figure; 
bode(Gyr1);
title('Bode Plot (k_p=0.01)');


%Case 2: kp=0.1:
% P controller 
kp  = 0.1; 
ki  = 0; %0.005;
C   = kp + ki/s;

L   = C*P2*P1;        % Loop transfer function 
Gyr2 = feedback(L,1);  % closed transfer function from r to y

% step response
figure; 
step(Gyr2);
title('Step Response (k_p=0.1)');


% freqency response 
figure; 
bode(Gyr2);
title('Bode Plot (k_p=0.1)');

% The case with kp = 0.01 results in a step response that reaches steady
% state without as many oscillations as the case with kp = 0.1. For the 
% frequency reponses, both cases have a phase decrease from 0 to -180
% degrees and overall, the shapes of both the magnitude and frequency
% responses are the same. However, when kp has a higher value, the
% magnitude response has more of a peak like shape at the crossover
% frequency before dipping and the phase response has a more obvious
% dip shape at this frequency. 

%% 4d)
% Dynamics 
s = tf('s');

P1 = b/(m*s+c);  % transfer function for tau to v
P2 = a*T/(s+a);  % transfer function for u to tau

% PI controller 
kp  = 0.1; 
ki  = 0.005;
C   = kp + ki/s;

L   = C*P2*P1;        % Loop transfer function 
Gyr3 = feedback(L,1);  % closed transfer function from r to y
disp(Gyr3);
% step response
figure; step(Gyr3)

% freqency response 
figure; bode(Gyr3)

%When integral control is added, the step response eventually reaches the
%steady state value instead of just oscillating around that region of the
%graph as it did in proportional control. At lower frequencies, the value
%of the magnitude response is larger (only a little larger in this case due
%to the small value of ki). Also, for both the magnitude and phase
%responses, the crossover frequency increases resulting in a wider
%magnitude and phase plot.

%% 4e)

%Case 1: kp = 0.01:
% ----------------------------------------------
%         Time domain simulation
% ----------------------------------------------

% Case 1: peice-wise constant reference
t = 0:0.01:500;
r = zeros(length(t),1);
for i = 1:length(t)   
    if t(i) <= 100
        r(i) = 10;
    elseif t(i) <= 200
        r(i) = 15;
    elseif t(i) <= 400
        r(i) = 13;
    else 
        r(i) = 10;
    end
end
y = lsim(Gyr1,r,t,0);
figure; plot(t,r,t,y);
title('Piecewise Reference Tracking (k_p=0.01)');
legend('Reference', 'Velocity');
set(gcf,'Position',[100 100 700 300])
    
% Case 2: periodic references
for period = [60, 50, 40, 20]   
    t = 0:0.01:6*period;        % time range -- 6 periods
    r = 10 + 2*sin(2*pi/period.*t);
    y = lsim(Gyr1,r,t,0);
    figure; plot(t,r,t,y);
    title('Periodic Reference Tracking with k_p=0.01 and Period =', num2str(period));
    set(gcf,'Position',[100 100 700 300])
    h = legend('Reference $r(t)$','Velocity $v(t)$',...
        'location','southeast','Interpreter','latex');
    set(h,'box','off');
    xlabel('Time $t$','Interpreter','latex');
    ylabel('Velocity $v$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex','fontsize',10);
end


%Case 2: kp = 0.1:
% ----------------------------------------------
%         Time domain simulation
% ----------------------------------------------

% Case 1: peice-wise constant reference
t = 0:0.01:500;
r = zeros(length(t),1);
for i = 1:length(t)   
    if t(i) <= 100
        r(i) = 10;
    elseif t(i) <= 200
        r(i) = 15;
    elseif t(i) <= 400
        r(i) = 13;
    else 
        r(i) = 10;
    end
end
y = lsim(Gyr2,r,t,0);
figure; plot(t,r,t,y);
title('Piecewise Reference Tracking (k_p=0.1)');
set(gcf,'Position',[100 100 700 300])
    
% Case 2: periodic references
for period = [60, 50, 40, 20]   
    t = 0:0.01:6*period;        % time range -- 6 periods
    r = 10 + 2*sin(2*pi/period.*t);
    y = lsim(Gyr2,r,t,0);
    figure; plot(t,r,t,y);
    
    set(gcf,'Position',[100 100 700 300])
    h = legend('Reference $r(t)$','Velocity $v(t)$',...
        'location','southeast','Interpreter','latex');
    set(h,'box','off');
    xlabel('Time $t$','Interpreter','latex');
    ylabel('Velocity $v$','Interpreter','latex');
    title('Periodic Reference Tracking with k_p=0.1 and Period =', num2str(period));
    set(gca,'TickLabelInterpreter','latex','fontsize',10);
end


%Case 3: PI Control
% ----------------------------------------------
%         Time domain simulation
% ----------------------------------------------

% Case 1: peice-wise constant reference
t = 0:0.01:500;
r = zeros(length(t),1);
for i = 1:length(t)   
    if t(i) <= 100
        r(i) = 10;
    elseif t(i) <= 200
        r(i) = 15;
    elseif t(i) <= 400
        r(i) = 13;
    else 
        r(i) = 10;
    end
end
y = lsim(Gyr3,r,t,0);
figure; plot(t,r,t,y);
title('Piecewise Reference Tracking (PI Control)');
set(gcf,'Position',[100 100 700 300])
    
% Case 2: periodic references
for period = [60, 50, 40, 20]   
    t = 0:0.01:6*period;        % time range -- 6 periods
    r = 10 + 2*sin(2*pi/period.*t);
    y = lsim(Gyr3,r,t,0);
    figure; plot(t,r,t,y);
    
    set(gcf,'Position',[100 100 700 300])
    h = legend('Reference $r(t)$','Velocity $v(t)$',...
        'location','southeast','Interpreter','latex');
    set(h,'box','off');
    xlabel('Time $t$','Interpreter','latex');
    ylabel('Velocity $v$','Interpreter','latex');
    title('Periodic Reference Tracking with PI Control and Period =', num2str(period));
    set(gca,'TickLabelInterpreter','latex','fontsize',10);
end


%Proportional control with a small kp does not accurately track the
%reference. The output velocity is much further away than the reference
%velocity. When the proportional constant is increased, the controller does
%a better job at tracking the reference with the velocity values being
%closer to the reference. However, this case isn't nearly as accurate the
%PI controller case. Adding integral control allows the controller to track
%the reference much more accurately and at some points, almost exactly.


%% 4f)
%To ensure safety in this system, it is important to have accurate
%reference tracking. For this, the kp and ki values must be properly
%adjusted. After experimenting with this simulation, it can be seen that
%proportional control is accurate for higher values of kp. However, one
%disadvantage of having higher values is the damping that will occur
%immediately after activating the system. In cruise control, this high
%damping will result in reference tracking errors that could lead to an
%accident. In order to get more accurate reference tracking with smaller
%proportional constants, it helps to have integral control. However, with
%integral control, higher constants will cause much more damping than lower
%constants, and it is much more sensitive to parameter changes than a
%proportional controller, so it's important to keep ki low. The best
%controller is one with a kp value high enough to avoid damping while also
%being able to perform reference tracking accurately and a ki value that is
%low enough to do the same. Finding this will require an optimization
%algorithm that finds the constants that result in the lowest amount of
%damping(possibly using peak difference).