
%ODE Simulations:

%Problem 4
%a)

%a
%--------------------------------------------------------------------------
%\dotx = 2x with x(0) = 1
%--------------------------------------------------------------------------
fa = @ (t, x)(2*x);
ode45(fa, [0, 10], 1);
xlabel('Time t');
ylabel('y');
title('ODE Simulation for equation (a)');

%--------------------------------------------------------------------------
% The graph of the solution displays exponential growth. As the t value 
% increases, the value of y increases exponentially. 
%--------------------------------------------------------------------------

%b
%--------------------------------------------------------------------------
%\dotx = -2x with x(0) = 1
%--------------------------------------------------------------------------
figure;
fb = @ (t, x)(-2*x);
ode45(fb, [0, 10], 1);
xlabel('Time t');
ylabel('y');
title('ODE Simulation for equation (b)');

%--------------------------------------------------------------------------
% The graph of the solution displays exponential decay. As the t value 
% increases, the value of y decreases exponentially. The y value at t = 0
% is 1 and from there it exponentially falls down until it approaches 0.
%-----------------------------------------------------------------------

%c
%--------------------------------------------------------------------------
%\dotx1 = x2 with x1(0) = 1
%\dotx2 = -x1-x2 with x2(0) = 1
%--------------------------------------------------------------------------
function dotx = fc(t, x)
    dotx1 = x(2);
    dotx2 = -x(1)-x(2);
    dotx = [dotx1; dotx2];
end

x = [1; 1];
figure;
ode45(@fc, [0, 10], x);
xlabel('Time t');
ylabel('y');
legend('x1', 'x2');
title('ODE Simulation for equation (c)');

%--------------------------------------------------------------------------
% The graph of the solution x1(t) shows x1(t) taking a value of 1 at t = 0
% and then momentarily rising up before slowly falling down below 0 and
% then adjusting itself back to 0 by time t = 10. The graph of the solution
% x2(t) shows x2(t) taking a value of 1 at t = 0 and steadily falling down
% to a negative value before it steadily rises up to slightly above 0. The
% graph then adjusts itself back to 0 by time t = 10. 
%--------------------------------------------------------------------------

%d
%--------------------------------------------------------------------------
%\dotx1 = x2 with x1(0) = 1
%\dotx2 = x1-x2 with x2(0) = 1
%--------------------------------------------------------------------------
function dotdx = fd(t, x)
    dotx1 = x(2);
    dotx2 = x(1)-x(2);
    dotdx = [dotx1; dotx2];
end

figure;
ode45(@fd, [0, 10], x);
xlabel('Time t');
ylabel('y');
legend('x1', 'x2');
title('ODE Simulation for equation (d)');

%--------------------------------------------------------------------------
% The graph of the solutions for x1(t) and x2(t) show exponential growth
% from 0. While both graphs have the same starting point, the graph for x1
% shows more exponential growth ending at a higher point than the graph for
% x2.
%--------------------------------------------------------------------------

%e
%--------------------------------------------------------------------------
%\dotx1 = x2 with x1(0) = 1
%\dotx2 = -x1 with x2(0) = 1
%--------------------------------------------------------------------------
function dotex = fe(t, x)
    dotx1 = x(2);
    dotx2 = -x(1);
    dotex = [dotx1; dotx2];
end

figure;
ode45(@fe, [0, 10], x);
xlabel('Time t');
ylabel('y');
legend('x1', 'x2');
title('ODE Simulation for equation (e)');

%--------------------------------------------------------------------------
% Both the graphs for x1(t) and x2(t) appear to be sinusoidal curves that
% are shifted versions of each other. 
%--------------------------------------------------------------------------


%d)
%--------------------------------------------------------------------------
%\dotx + 3x = 5 with x(0) = 2
%--------------------------------------------------------------------------
figure;
f4d = @(t, x)(5-3*x);
ode45(f4d, [0, 10], 2);
xlabel('Time t');
ylabel('y');
title('ODE Simulation for problem(d)');

%Problem 5
%a)
%--------------------------------------------------------------------------
%\dotx = sigma(y-x)
%\doty = x(rho-z)-y
%\dotz = xy-beta*z
%--------------------------------------------------------------------------
function dot5x = f5a(t, x)
    sig = 10;
    rho = 28;
    beta = 8 / 3;
    dotx1 = sig * (x(2) - x(1));
    dotx2 = x(1) * (rho - x(3)) - x(2);
    dotx3 = (x(1) * x(2)) - (beta * x(3));
    dot5x = [dotx1; dotx2; dotx3];
end

%1)
x01 = [1; 1; 1];
figure;
ode45(@f5a, [0, 30], x01);
xlabel('Time t');
ylabel('y');
title('ODE Simulation with IC [x0,y0,z0] = [1,1,1]');
[ts51, ys51] = ode45(@f5a, [0, 30], x01);


%2)
x02 =  [1.01; 1; 1];
figure;
ode45(@f5a, [0, 30], x02);
xlabel('Time t');
ylabel('y');
title('ODE Simulation with IC [x0,y0,z0] = [1.01,1,1]');
[ts52, ys52] = ode45(@f5a, [0, 30], x02);


%b)
figure;
plot3(ys51(:,1), ys51(:,2), ys51(:,3), 'b', 'LineWidth', 1.5)
hold on
plot3(ys52(:,1), ys52(:,2), ys52(:,3), 'r', 'LineWidth', 1.5)
legend('IC with all 1s', 'IC with x0=1.01');
xlabel('x(t)');
ylabel('y(t)');
zlabel('z(t)');
title('Lorenz Attractor Simulation');

%--------------------------------------------------------------------------
% The figure that is shown in the 3d graph is an attractor, which is a
% fractal space that consists of butterfly shaped curves representing the 
% outputs of the differential equations with the different initial
% conditions. Although both the initial conditions may be closely spaced,
% there is a high level of divergence in the paths of both graphs and although they
% are bounded to certain values in the x, y, and z axes, the paths are
% unpredictable with a spiraling structure.
%--------------------------------------------------------------------------

%Problem 6
%a)
%open loop
f_resp_open = @(beta) 1 ./ beta;
beta = linspace(0, 30, 100);
v_resp_open = f_resp_open(beta);
figure;
plot(beta, v_resp_open);
title('Open Loop Steady State');
xlabel('beta');
ylabel('v resp Open Loop');

%feedback
f_resp_fb = @(beta) 10 ./ (beta + 10);
v_resp_fb = f_resp_fb(beta);
figure;
plot(beta, v_resp_fb);
title('Feedback Steady State');
xlabel('beta');
ylabel('v resp Feedback');
