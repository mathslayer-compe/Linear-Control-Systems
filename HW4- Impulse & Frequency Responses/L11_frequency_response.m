%% 3a)

% -------------------------------------------------------------------------
%
%  ECE 171A: Linear Control System Theory
%  Frequency response - time plot
%
% -------------------------------------------------------------------------


close all; clear;

% --------------------------------------------------
%    state-space system
% --------------------------------------------------
k = 1;
m = 1;
c = 0.2;
A1   = [0 1; -k/m, -c/m]; 
B1   = [0; 1/m]; 
C1   = [1, 0];
D1   = 0;
sys1 = ss(A1, B1, C1, D1);  % stable system 1
x0 = [0; 0];

% --------------------------------------------------
% Case 1: input signal (sin(0.1t))   
% --------------------------------------------------

omega = 0.1;
T     = 2*pi/omega*10; 
t     = 0:.0001:T;

u1    = sin(omega*t);
y1    = lsim(sys1,u1,t,x0);

figure;  FontSize = 10;
plot(t,u1,'black',t,y1,'b');
h = legend('Input','output for system','Interpreter','latex');
set(h,'box','off')
title('Frequency response for $\sin(0.1t)$','Interpreter','latex','fontsize',FontSize);
axis([0,T, -3, 3]);
grid on;
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);
set(gcf,'Position',[100 100 500 300])
print(gcf,'L11_fre2','-painters','-depsc','-r300')
[pks_y, locs_y] = findpeaks(y1, t);
[pks_u, locs_u] = findpeaks(u1, t);
Ay = mean(pks_y);
Au = mean(pks_u);
M1 = Ay / Au;
n = min(length(locs_y), length(locs_u));
T_vector = locs_y(1:n) - locs_u(1:n);
deltaT = mean(T_vector);
period = 2*pi*omega;
phase1 = (-2*pi*deltaT) / (T/10);

% --------------------------------------------------
% Case 2: input signal (sin(0.5t))   
% --------------------------------------------------

omega = 0.5;
T     = 2*pi/omega*10; 
t     = 0:.0001:T;

u1    = sin(omega*t);
y1    = lsim(sys1,u1,t, x0);

figure;  FontSize = 10;  
plot(t,u1,'black',t,y1,'b');
h = legend('Input','output for system','Interpreter','latex');
set(h,'box','off')
title('Frequency response for $\sin(0.5t)$','Interpreter','latex','fontsize',FontSize);
axis([0,T, -3, 3]);
grid on;
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);
set(gcf,'Position',[100 100 500 300])
print(gcf,'L11_fre2','-painters','-depsc','-r300')
[pks_y, locs_y] = findpeaks(y1, t);
[pks_u, locs_u] = findpeaks(u1, t);
Ay = mean(pks_y);
Au = mean(pks_u);
M2 = Ay / Au;
n = min(length(locs_y), length(locs_u));
T_vector = locs_y(1:n) - locs_u(1:n);
deltaT = mean(T_vector);
period = 2*pi*omega;
phase2 = (-2*pi*deltaT) / (T/10);


% --------------------------------------------------
% Case 3: input signal (sin(t))   
% --------------------------------------------------

omega = 1;
T     = 2*pi/omega*10; 
t     = 0:.0001:T;

u1    = sin(omega*t);
y1    = lsim(sys1,u1,t, x0);

figure;  FontSize = 10;  
plot(t,u1,'black',t,y1,'b');
h = legend('Input','output for system','Interpreter','latex');
set(h,'box','off')
title('Frequency response for $\sin(t)$','Interpreter','latex','fontsize',FontSize);
axis([0,T, -3, 3]);
grid on;
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);
set(gcf,'Position',[100 100 500 300])
print(gcf,'L11_fre2','-painters','-depsc','-r300')
[pks_y, locs_y] = findpeaks(y1, t);
[pks_u, locs_u] = findpeaks(u1, t);
Ay = mean(pks_y);
Au = mean(pks_u);
M3 = Ay / Au;
n = min(length(locs_y), length(locs_u));
T_vector = locs_y(1:n) - locs_u(1:n);
deltaT = mean(T_vector);
period = 2*pi*omega;
phase3 = (-2*pi*deltaT) / (T/10);

% --------------------------------------------------
% Case 4: input signal (sin(1.5t))   
% --------------------------------------------------

omega = 1.5;
T     = 2*pi/omega*10; 
t     = 0:.0001:T;

u1    = sin(omega*t);
y1    = lsim(sys1,u1,t, x0);

figure;  FontSize = 10;  
plot(t,u1,'black',t,y1,'b');
h = legend('Input','output for system','Interpreter','latex');
set(h,'box','off')
title('Frequency response for $\sin(1.5t)$','Interpreter','latex','fontsize',FontSize);
axis([0,T, -3, 3]);
grid on;
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);
set(gcf,'Position',[100 100 500 300])
print(gcf,'L11_fre2','-painters','-depsc','-r300')
[pks_y, locs_y] = findpeaks(y1, t);
[pks_u, locs_u] = findpeaks(u1, t);
Ay = mean(pks_y);
Au = mean(pks_u);
M4 = Ay / Au;
n = min(length(locs_y), length(locs_u));
T_vector = locs_y(1:n) - locs_u(1:n);
deltaT = mean(T_vector);
period = 2*pi*omega;
phase4 = (-2*pi*deltaT) / (T/10);

% --------------------------------------------------
% Case 5: input signal (sin(2t))   
% --------------------------------------------------

omega = 2;
T     = 2*pi/omega*10; 
t     = 0:.0001:T;

u1    = sin(omega*t);
y1    = lsim(sys1,u1,t, x0);

figure;  FontSize = 10;  
plot(t,u1,'black',t,y1,'b');
h = legend('Input','output for system','Interpreter','latex');
set(h,'box','off')
title('Frequency response for $\sin(2t)$','Interpreter','latex','fontsize',FontSize);
axis([0,T, -3, 3]);
grid on;
set(gca,'TickLabelInterpreter','latex','fontsize',FontSize);
set(gcf,'Position',[100 100 500 300])
print(gcf,'L11_fre2','-painters','-depsc','-r300')
[pks_y, locs_y] = findpeaks(y1, t);
[pks_u, locs_u] = findpeaks(u1, t);
Ay = mean(pks_y);
Au = mean(pks_u);
M5 = Ay / Au;
n = min(length(locs_y), length(locs_u));
T_vector = locs_y(1:n) - locs_u(1:n);
deltaT = mean(T_vector);
period = 2*pi*omega;
phase5 = (-2*pi*deltaT) / (T/10);

%% 3b)
mag = [M1; M2; M3; M4; M5];
phi = [phase1; phase2; phase3; phase4; phase5];
omegas = [0.1; 0.5; 1; 1.5; 2];
tabl = table(omegas, 20*log10(mag), rad2deg(phi), 'VariableNames', {'omega', 'Gain M', 'Phase phi'});
disp(tabl);

%% 3c)
figure;
loglog(omegas, 20*log10(mag));
title('Magnitude M as a function of Frequency on log-log scale');
xlabel('Frequency (rad/s)');
ylabel('Magnitude M (dB)');

figure;
semilogx(omegas, rad2deg(phi));
title('Phase as a function of Frequency on log-linear scale');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');

%% 4c) 
omega0 = 1;
%--------------------------------------------------------------------------
% Case 1: damping = 0.1
%--------------------------------------------------------------------------
damping = 0.1;
k1 = (4 * damping * omega0) - 8;
k2 = (2*omega0^2) - (4 * damping * omega0) + 6;
A = [-3 2;1, -1];
B = [0.5;0];
C = [0 1];
D = 0;
K1 = [k1 k2];
kr1 = 2*omega0^2;
sys = ss(A-B.*K1, kr1*B, C, D);
e1 = eig(A-B.*K1);
disp('Eigenvalues (damping = 0.1):');
disp(e1);
fprintf('kr(damping = 0.1): = %d\n', kr1);
disp('K matrix (damping=0.1): ')
disp(K1);
figure;
step(sys);

%--------------------------------------------------------------------------
% Case 2: damping = 0.4
%--------------------------------------------------------------------------
damping = 0.4;
k1 = (4 * damping * omega0) - 8;
k2 = (2*omega0^2) - (4 * damping * omega0) + 6;
A = [-3 2;1, -1];
B = [0.5;0];
C = [0 1];
D = 0;
K2 = [k1 k2];
kr2 = 2*omega0^2;
sys = ss(A-B.*K2, kr2*B, C, D);
e2 = eig(A-B.*K2);
disp('Eigenvalues (damping = 0.4):');
disp(e2);
fprintf('kr(damping = 0.4): = %d\n', kr2);
disp('K matrix (damping=0.4): ')
disp(K2);
figure;
step(sys);

%--------------------------------------------------------------------------
% Case 3: damping = 0.7
%--------------------------------------------------------------------------
damping = 0.7;
k1 = (4 * damping * omega0) - 8;
k2 = (2*omega0^2) - (4 * damping * omega0) + 6;
A = [-3 2;1, -1];
B = [0.5;0];
C = [0 1];
D = 0;
K3 = [k1 k2];
kr3 = 2*omega0^2;
sys = ss(A-B.*K3, kr3*B, C, D);
e3 = eig(A-B.*K3);
disp('Eigenvalues (damping = 0.7):');
disp(e3);
fprintf('kr(damping = 0.7): = %d\n', kr3);
disp('K matrix (damping=0.7): ')
disp(K3);
figure;
step(sys);

%--------------------------------------------------------------------------
% Case 4: damping = 0.9
%--------------------------------------------------------------------------
damping = 0.9;
k1 = (4 * damping * omega0) - 8;
k2 = (2*omega0^2) - (4 * damping * omega0) + 6;
A = [-3 2;1, -1];
B = [0.5;0];
C = [0 1];
D = 0;
K4 = [k1 k2];
kr4 = 2*omega0^2;
sys = ss(A-B.*K4, kr4*B, C, D);
e4 = eig(A-B.*K4);
disp('Eigenvalues (damping = 0.9):');
disp(e4);
fprintf('kr(damping = 0.9): = %d\n', kr4);
disp('K matrix (damping=0.9): ')
disp(K4);
figure;
step(sys);

%--------------------------------------------------------------------------
% Observation:
% The amount of oscillations that the system has to go through before
% reaching steady state decreases as the damping constant increases from
% 0.1 to 0.9.
%--------------------------------------------------------------------------

%% 4d)
%--------------------------------------------------------------------------
% I would advise healthcare providers to develop and dose controllers whose
% transient overshoots don't exceed the upper limit of healthy nutrient
% intake in order to not induce any toxicity in the bloodstream and to
% develop pricing models based on individual patient data that show the
% dosage that results in the greatest health gain per dollar spent.
% Socially, there are people with different income levels, but all of them
% should get access to nutrition education so that they can learn how to
% live healthy lifestyles from a young age, and there should be publically
% funded healthcare programs that allow patients to get access to
% healthcare for a minimal cost or maybe even for free. To address the
% problem of each individual having different patterns, the solution is the
% pricing model as mentioned before that is thoroughly simulated in
% multiple scenarios with the patient's parameters in order to find the
% optimal price/dosage.
%--------------------------------------------------------------------------