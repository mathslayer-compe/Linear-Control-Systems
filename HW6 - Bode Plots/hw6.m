%% 1d)
k = 20;
C = k;
den = [1 8 9 -18];
P = tf(1, den);
Gyr = feedback(P*C, 1);
disp(Gyr);
figure;
pzmap(Gyr);

%The poles of this system are all on the left half plane, which indicates
%that the system is stable. According to my answer in part c, in order for
%a system to be stable, the range of k must be 18<k<90. Since k=20 falls
%within the range, the system is expected to be stable.

figure;
step(Gyr);

%% 2f)
a = 1;
b = 1;
omega = 1;
zeta = 0.5;
tau = 1;

%% 2a)
num = [b a*b];
den = [1 0];
sys = tf(num, den);
figure;
bode(sys);

%% 2b)
num = [1 a];
den = [1 100*a];
sys = tf(num, den);
figure;
bode(sys);

%% 2c)
num = [1 100*a];
den = [1 a];
sys = tf(num, den);
figure;
bode(sys);

%% 2d)
num = 1;
den = [1 2*zeta*omega omega^2 0];
sys = tf(num, den);
figure;
bode(sys);

%% 2e)
sys = tf(1, 1, 'InputDelay', tau);
figure;
bode(sys);

% These plots are consistent with my sketches.

%% 3c)
num = 1;
den = [1 3 2];
sys = tf(num, den);
figure;
margin(sys);

%The gain margin is the maximum gain increase or decrease of a system that
%doesn't compromise stability. This system is already stable because it
%has poles with negative real parts, so the gain margin will be infinity.
%This same transfer function block was part of the open loop system in
%problem 1b except that it was multiplied by a constant k. The range of the
%constant k was determined to be k>-2 for stability, and this aligns with
%the infinite gain margin of L(s) because it indicates that we can increase
%k from -2 to any value without losing stability.


