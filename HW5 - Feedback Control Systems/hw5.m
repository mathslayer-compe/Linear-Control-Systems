%% 1d) 
A = [-1 -2; 1, 0];
B = [1; 0];
C = [0 1];
D = 0;
sys = ss(A, B, C, D);
figure;
step(sys);

%Verification of Analytical Computation:
num = 1;
den = [1 1 2];
sys = tf(num, den);
figure;
step(sys);

%Since both the step responses are the same, the analytical computation of
%the transfer function is consistent with the result of the numerical
%simulation.

%% 2c)
m = 1;
k = 1;
c = 0.2;
A = [0 1; -k/m, -c/m];
B = [0; 1/m];
C = [1 0];
D = 0;
sys = ss(A, B, C, D);
figure;
bode(sys, {0.01, 10});

%The numbers in my table, although in different units, are consistent with
%the Bode plot numbers.

%% 4b)
m1 = 1;
c1 = 1;
k1 = 1;
m2 = 1;
k2 = 1;

num1 = [m2 0 k2];
den = [m1*m2 m2*c1 ((m1*k2)+m2*(k1+k2)) k2*c1 k2*(k1+k2)-k2^2];
num2 = k2;
sys1 = tf(num1, den);
figure;
pzmap(sys1);
sys2 = tf(num2, den);
figure;
pzmap(sys2);

fprintf('Poles q_1:\n')
disp(pole(sys1));
fprintf('Zeros q_1:\n')
disp(zero(sys1));
fprintf('Poles q_2:\n')
disp(pole(sys2));
fprintf('Zeros q_1:\n')
disp(zero(sys2));

%The poles in the pole zero map correspond to the location of the points
%marked with an x while the zeros correspond to the location of the points
%marked with o. The poles of both the systems are the same as they share
%the same denominator, and they are all to the left of Real = 0.

%% 4c) 
figure;
bode(sys1);

T  = 100;        
dt = 0.01;
t  = 0:dt:T;
w = [0.578, 1.0, 1.1];
u1 = sin(w(1)*t);
u2 = sin(w(2)*t);
u3 = sin(w(3)*t);
x0 = [0 0 0 0]; 
y1 = lsim(sys1, u1, t, x0);
y2 = lsim(sys1, u2, t, x0);
y3 = lsim(sys1, u3, t, x0);
figure;
plot(t,y1,'b', t,y2,'r', t,y3,'k')
legend('\omega=0.578','\omega=1.0','\omega=1.1')
xlabel('Time (s)')
ylabel('q_1(t)')
title('Response of q_1 to sin(\omega t)')

%These responses are consistent with the Bode plots. The steady state
%amplitudes line up with the Bode magnitude curve. The numerator of this
%system is essentially s^2+1, which means that the system has zeros at
%s=+/-i\omega with \omega=1. At \omega=1, a sinsuoidal input will not
%produce a long term output in this system and will go down to 0. In a
%physical system, this can be seen as a tuned damper canceling out
%the motion of a primary mass at the natural frequency.
