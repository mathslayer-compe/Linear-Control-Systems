%% 1d
% ------------------------------------------------------------------------ 
%
%  ECE171A: HW3 Problem 1 - sample code
%  Original Writen by Yang Zheng
%  Spring, 2025
%  Simulation of Parameters by Akhil Nallacheruvu
%
% ------------------------------------------------------------------------ 
close all;

% parameters
r = 0.1;
e = 0.1;
c = 100;
b = 0.2; 
a = 0.5;
w = -0.61;


% dynamics
f = @(t,x) [r*x(1) - a*x(2); b*x(1)*x(2)/(c+x(1))-e*x(2)]; 

H = linspace(0,200,40);
G = linspace(0,40,20);

[x,y] = meshgrid(H,G);
u = zeros(size(x));
v = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0;   % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

% We use the quiver command to plot our vector field
figure; quiver(x,y,u,v,'r'); 
xlabel('$H$: Hare','Interpreter','latex')
ylabel('$G$: Tiger','Interpreter','latex')
axis tight equal;

% Plot some trajectories
hold on

H0 = 10:10:200;

for i=1:length(H0)
    [ts,ys] = ode45(f,[0,25],[H0(i);0.2*H0(i)]);     % ode45 simulations
    plot(ys(:,1),ys(:,2),'linewidth',1.2)
    plot(ys(1,1),ys(1,2),'bo','MarkerSize',10,'LineWidth',1.2)        % starting point
    plot(ys(end,1),ys(end,2),'ks','MarkerSize',10,'LineWidth',1.2)    % ending point
end

ylim([0,40]); xlim([0,200]);
hold off
set(gcf,'Position',[150 150 900 250])

%% 1f

% dynamics
ff = @(t,x) [r*x(1) - a*x(2) + w*(x(1)-100); b*x(1)*x(2)/(c+x(1))-e*x(2)];

% We use the quiver command to plot our vector field
figure; quiver(x,y,u,v,'r'); 
xlabel('$H$: Hare','Interpreter','latex')
ylabel('$G$: Tiger','Interpreter','latex')
axis tight equal;

% Plot some trajectories
hold on

for i=1:length(H0)
    [ts,ys] = ode45(ff,[0,500],[H0(i);0.2*H0(i)]);     % ode45 simulations
    plot(ys(:,1),ys(:,2),'linewidth',1.2)
    plot(ys(1,1),ys(1,2),'bo','MarkerSize',10,'LineWidth',1.2)        % starting point
    plot(ys(end,1),ys(end,2),'ks','MarkerSize',10,'LineWidth',1.2)    % ending point
end

ylim([0,40]); xlim([0,200]);
hold off
set(gcf,'Position',[150 150 900 250])

% Comparison:
% The phase portrait in part (d) shows trajectories from each initial point
% moving away from the equilibrium point. The phase portrait looks like a 
% source. This happens because the system is unstable at the equilibrium
% point without any Jacobian linearization. The phase portrait in part (f)
% is much closer to a sink with the trajectories from the initial points
% spinning more directly near to the equilibrium point. This happens
% because the linearized system is stable at the equilibrium point.

%% 2a)
f2a_u1 = @(t,x) [x(2); sin((2*pi)/4)-(0.2*x(2))-x(1)];
figure;
ode45(f2a_u1,[0,100],[0;0]);
xlabel('t');
ylabel('x');
legend('$x_1$', '$x_2$', 'Interpreter', 'Latex');
title('Spring Mass Dynamics for Input $u_1$', 'Interpreter', 'Latex');


f2a_u2 = @(t,x) [x(2); sin((2*pi)/20)-(0.2*x(2))-x(1)];
figure;
ode45(f2a_u2, [0,100], [0;0]);
xlabel('t');
ylabel('x');
legend('$x_1$', '$x_2$', 'Interpreter', 'Latex');
title('Spring Mass Dynamics for Input $u_2$', 'Interpreter', 'Latex');

f2a_u3 = @(t,x) [x(2); sin((2*pi)/4)+sin((2*pi)/20)-(0.2*x(2))-x(1)];
figure;
ode45(f2a_u3, [0,100], [0;0]);
xlabel('t');
ylabel('x');
legend('$x_1$', '$x_2$', 'Interpreter', 'Latex');
title('Spring Mass Dynamics for Input $u_3$', 'Interpreter', 'Latex');

%% 2b)
%Observations:
%The points on the plot for the outputs of the 3rd input u3 are the sum of
%the points on the plot for the outputs of the 1st input u1 and the second
%input u2. This happens because the system is linear, which means that
%if an input is the sum of inputs, then the output will be the sum of the
%outputs of each individual input.