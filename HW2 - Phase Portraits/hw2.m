% Phase Portrait Plotting

%% 1b) 
% Define the variables
L = 0.1;
Vs = 1;
C = 0.2;
R = 1;

% Define the functions
f = @(t,x) [(-1/L) * x(2) + (1/L) * Vs; 
            (1/C) * x(1) - (1/(R*C)) * x(2)];

% Step 2: Create a grid of, e.g., 30x30 points.
y1 = linspace(-2,8,30);
y2 = linspace(-2,2,30);

% Step 3: creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. 
% Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 columns
[x,y] = meshgrid(y1,y2);

% Step 4: computing the vector field
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

% Step 5: we use the quiver command to plot our vector field

figure; quiver(x,y,u,v,'r'); 
xlabel('$x$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
axis tight equal;
set(gcf,'Position',[150 150 600 300])

% Step 6: Plotting solutions on the vector field
% Let's plot a few solutions on the vector field. 
% We will consider the solutions where y1(0)=0, and values of y2(0) = [0 0.5 1 1.5 2.1 2.5], 
% in otherwords, we start the pendulum at an angle of zero, with some angular velocity.
hold on
for y20 = [0 0.5 1 1.5 2.1 2.5] 
    [ts,ys] = ode45(f,[0,50],[0;y20]);     % ode45 simulations
    plot(ys(:,1),ys(:,2),'linewidth',1.2)
    plot(ys(1,1),ys(1,2),'bo','MarkerSize',10,'LineWidth',1.2)        % starting point
    plot(ys(end,1),ys(end,2),'ks','MarkerSize',10,'LineWidth',1.2)    % ending point
end
ylim([-2.5,2.5]); xlim([-2,8]);
hold off

%% 1c) 
fc = @(t,x) [( (t < 1)*0 + (t >= 1)*1 - x(2) )/L; 
            (1/C) * x(1) - (1/(R*C)) * x(2)];

figure;
ode45(fc, [0 10], [0 0]);
xlabel('t');
ylabel('x');
title('RLC Circuit State Trajectory in Time Domain');
legend('$\dot{x_1}$', '$\dot{x_2}$','Interpreter','latex');


%% 4a)
f4a = @(t,x)[x(2); -sin(x(1))];
u4a = zeros(size(x));
v4a = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t = 0;
for i = 1:numel(x)
    Yprime = f4a(t,[x(i); y(i)]);
    u4a(i) = Yprime(1);
    v4a(i) = Yprime(2);
end

% Step 5: we use the quiver command to plot our vector field

figure; quiver(x,y,u4a,v4a,'r'); 
xlabel('$y_1=\theta$','Interpreter','latex')
ylabel('$y_2=\dot{\theta}$','Interpreter','latex')
axis tight equal;
set(gcf,'Position',[150 150 600 300])


% Step 6: Plotting solutions on the vector field
% Let's plot a few solutions on the vector field. 
% We will consider the solutions where y1(0)=0, and values of y2(0) = [0 0.5 1 1.5 2.1 2.5], 
% in otherwords, we start the pendulum at an angle of zero, with some angular velocity.
hold on
for y20 = [0 0.5 1 1.5 2.1 2.5]
    [ts,ys] = ode45(f4a,[0,50],[0;y20]);     % ode45 simulations
    plot(ys(:,1),ys(:,2),'linewidth',1.2)
    plot(ys(1,1),ys(1,2),'bo','MarkerSize',10,'LineWidth',1.2)        % starting point
    plot(ys(end,1),ys(end,2),'ks','MarkerSize',10,'LineWidth',1.2)    % ending point
end
ylim([-2.5,2.5]); xlim([-2,8]);
hold off

%Figure Description:
% This figure displays a center around an equilibrium point when the 
% initial y points are 0, 0.5, 1, and 1.5. When the initial points are 2.1
% and 2.5, then the shape of the plot is a curve that doesn't approach the
% equilibrium point at any point in time.

%% 4b)
f4b = @(t,x)[x(2); (-0.2*x(2))-sin(x(1))];
u4b = zeros(size(x));
v4b = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t = 0;
for i = 1:numel(x)
    Yprime = f4b(t,[x(i); y(i)]);
    u4b(i) = Yprime(1);
    v4b(i) = Yprime(2);
end

% Step 5: we use the quiver command to plot our vector field

figure; quiver(x,y,u4b,v4b,'r'); 
xlabel('$y_1=\theta$','Interpreter','latex')
ylabel('$y_2=\dot{\theta}$','Interpreter','latex')
axis tight equal;
set(gcf,'Position',[150 150 600 300])

% Step 6: Plotting solutions on the vector field
% Let's plot a few solutions on the vector field. 
% We will consider the solutions where y1(0)=0, and values of y2(0) = [0 0.5 1 1.5 2.1 2.5], 
% in otherwords, we start the pendulum at an angle of zero, with some angular velocity.
hold on
for y20 = [0 0.5 1 1.5 2.1 2.5]
    [ts,ys] = ode45(f4b,[0,50],[0;y20]);     % ode45 simulations
    plot(ys(:,1),ys(:,2),'linewidth',1.2)
    plot(ys(1,1),ys(1,2),'bo','MarkerSize',10,'LineWidth',1.2)        % starting point
    plot(ys(end,1),ys(end,2),'ks','MarkerSize',10,'LineWidth',1.2)    % ending point
end
ylim([-2.5,2.5]); xlim([-2,8]);
hold off

%Figure Description:
% This figure displays two sinks. This indicates that both the equilibrium
% points in this system are asymptotically stable.

%% 4c)
f4c = @(t,x)[x(2); (-0.5*x(2))-sin(x(1))];
u4c = zeros(size(x));
v4c = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t = 0;
for i = 1:numel(x)
    Yprime = f4c(t,[x(i); y(i)]);
    u4c(i) = Yprime(1);
    v4c(i) = Yprime(2);
end

% Step 5: we use the quiver command to plot our vector field

figure; quiver(x,y,u4c,v4c,'r'); 
xlabel('$y_1=\theta$','Interpreter','latex')
ylabel('$y_2=\dot{\theta}$','Interpreter','latex')
axis tight equal;
set(gcf,'Position',[150 150 600 300])

% Step 6: Plotting solutions on the vector field
% Let's plot a few solutions on the vector field. 
% We will consider the solutions where y1(0)=0, and values of y2(0) = [0 0.5 1 1.5 2.1 2.5], 
% in otherwords, we start the pendulum at an angle of zero, with some angular velocity.
hold on
for y20 = [0 0.5 1 1.5 2.1 2.5]
    [ts,ys] = ode45(f4c,[0,50],[0;y20]);     % ode45 simulations
    plot(ys(:,1),ys(:,2),'linewidth',1.2)
    plot(ys(1,1),ys(1,2),'bo','MarkerSize',10,'LineWidth',1.2)        % starting point
    plot(ys(end,1),ys(end,2),'ks','MarkerSize',10,'LineWidth',1.2)    % ending point
end
ylim([-2.5,2.5]); xlim([-2,8]);
hold off

%Figure Description:
% This figure displays a sink. This indicates that the equilibrium point is
% asymptotically stable.