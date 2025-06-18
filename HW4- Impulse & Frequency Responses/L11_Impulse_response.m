%% 2b)
% -------------------------------------------------------------------------
%
%  ECE 171A: Linear Control System Theory
%  Impulse response - time plot
%
% -------------------------------------------------------------------------
clc;
close all; clear; 

% --------------------------------------------------
% Case 1: u(t) = p_e(t)
% --------------------------------------------------
% --------------------------------------------------
%    state-space system
% --------------------------------------------------
A    = [-1 1; 0, -2]; 
B    = [0; 1]; 
C    = [0, 1];
D    = 0;
sys1 = ss(A, B, C, D);   % A stable system

% --------------------------------------------------
% define the impulse input signal 
%   u and simulate the impuse response
% --------------------------------------------------

% define the impulse signal (which is an approximation)
T   = 7;
t   = 0:.0001:T;
eps = 0.01;
u1  = zeros(size(t));
i   = 1;
while t(i) < eps
    u1(i)= 1/eps;
    i    = i+1;
end

% --------------------------------------------------
% simulate impulse response
% --------------------------------------------------
x0   = [0;0];              % zero initial condition
y1   = lsim(sys1,u1,t,x0); % stable system

% --------------------------------------------------
% plot the input signal and output signal
% --------------------------------------------------
figure; FontSize = 8;
subplot(2,1,1); plot(t,u1);    % input signal
axis([-0.2,T 0 1/eps*2]);
xlabel('time t','Interpreter','latex');
ylabel('input signal','Interpreter','latex');
title('Impulse signal $u(t) = p_{\epsilon}(t)$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

subplot(2,1,2);           % output signal
plot(t,y1);
axis([-0.25, T, 0, 1]);
title('Stable system','Interpreter','latex');
xlabel('time t','Interpreter','latex');
ylabel('Output signal $y$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);


set(gcf,'Position',[100 100 500 400])

% --------------------------------------------------
% Case 2: u(t) = p_e(t-1)
% --------------------------------------------------
tau  = 1;
i    = 1;
u3   = zeros(size(u1));
u3(tau/0.0001:end) = u1(1:end-tau/0.0001+1); % a shift of the input

y3   = lsim(sys1,u3,t,x0);

figure;
subplot(2,1,1); plot(t,u3); % input signal
axis([-0.2 T 0 1/eps*2]);
xlabel('Time $t$','Interpreter','latex');
ylabel('Input signal','Interpreter','latex');
title('Impulse signal $u(t) = p_{\epsilon}(t-1)$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

subplot(2,1,2); % output signal
plot(t,y3);
axis([-0.25 T 0 1]);
title('Impulse response under $u(t) = p_{\epsilon}(t-1)$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

set(gcf,'Position',[100 100 500 400])



% --------------------------------------------------
% Case 3: u(t) p_e(t-2)
% --------------------------------------------------

tau  = 2;
i    = 1;
u4   = zeros(size(u1));
u4(tau/0.0001:end) = u1(1:end-tau/0.0001+1); % a shift of the input

y4   = lsim(sys1,u4,t,x0);

figure;
subplot(2,1,1);  plot(t,u4); % input signal
axis([-0.2 T 0 1/eps*2]);
xlabel('Time $t$','Interpreter','latex');
ylabel('Input signal','Interpreter','latex');
title('Impulse signal $u(t) = p_{\epsilon}(t-2)$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

subplot(2,1,2); % output signal
plot(t,y4);
axis([-0.25 T 0 1]);
title('Impulse response under $u(t) = p_{\epsilon}(t-2)$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

set(gcf,'Position',[100 100 500 400])

% --------------------------------------------------
% Case 4: u(t) = p_e(t) + p_e(t-1) + p_e(t-2)
% --------------------------------------------------

% recall that u3 and u4 are shiftted signals of u1
u5 = u1 + u3 + u4;
y5 = lsim(sys1,u5,t,x0);

figure;
subplot(3,1,1);  plot(t,u5);   % input signal
axis([-0.2,T 0 1/eps*2]);
xlabel('Time $t$','Interpreter','latex');
ylabel('Input signal','Interpreter','latex');
title('Impulse signal $u(t) = p_{\epsilon}(t) + p_{\epsilon}(t-1) + p_{\epsilon}(t-2)$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);


subplot(3,1,2); plot(t,y5);% output signal
axis([-0.25,T 0 1.5]);
title('Impulse response under $u(t)$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

% let's compute y1+y2+y3 and compare it with y4
subplot(3,1,3);
plot(t,y1+y3+y4);
axis([-0.25,T 0 1.5]);
title('$y_1+y_3+y_4$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

set(gcf,'Position',[100 100 500 400])

%--------------------------------------------------------------------------
% Observations:
% In case 1, the input surges at time t=0 to 100 and almost immediately 
% falls back down to 0. The output exhibhits the same behavior at time t=0
% and surges up to the value y=1 and then exponentially decays to 0
% indicating that the system is stable. In case 2, the input has similar
% behavior to that in case 1 except it is delayed by 1 unit of time, so the
% input's surge from 0 to 100 happens at t=1 instead of t=0. This delay can
% also be seen in the output, which has similar behavior to the output in
% case 1. There is a surge to y=1 at time t=1 and then an exponential decay
% back down to. Case 3 has the same behavior except with a delay of 2 units
% of time. The surges and decreases/decays of the input and output occur at
% and after t=2. All 3 systems are stable. Case 4's input is the sum of the
% inputs in the previous cases, which can be seen on the input graph as 3
% signals (pulses) at times t=0,1, and 2. What's more interesting however,
% is that the output is also the sum of the outputs of the 3 previous
% cases. There is a pulse at time t=0 to y=1 and then an exponential decay
% until another pulse at time t=1 to y=1 and then an exponential decay and
% finally, there is a pulse at time t=2 and a continuous(uninterrupted) 
% exponential decay down to 0. This indicates that the system has linearity
% as the sum of individual inputs results in an output that is the sum of
% the outputs of the individual inputs.
%--------------------------------------------------------------------------

%% 2c)
u = zeros(size(t));
x0 = [0; 1];
y = lsim(sys1, u, t, x0);

figure; FontSize = 8;
subplot(2,1,1); plot(t,u);    % input signal
axis([-0.2,T 0 1/eps*2]);
xlabel('time t','Interpreter','latex');
ylabel('input signal','Interpreter','latex');
title('Impulse signal $u(t)=0$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

subplot(2,1,2);           % output signal
plot(t,y);
axis([-0.25, T, 0, 1]);
title('Stable system','Interpreter','latex');
xlabel('time t','Interpreter','latex');
ylabel('Output signal $y$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);

%--------------------------------------------------------------------------
% Comparison:
% There is no pulse in the input signal unlike case 1. It is consistently
% 0. The output signal starts with the value 1 and exponentially decays to 
% 0, and there is no pulse in the output signal either. However, just like
% case 1, the system is stable, so it has the exponential decay in common.
%--------------------------------------------------------------------------