%% 1d)
a = 1;
tau = 0.25;
num = a;
den = [tau a*tau 0];
P = tf(num, den);
P_del = tf(num, den, 'InputDelay', tau);
P = P - P_del;
figure;
nyquist(P);
numC = [6 15];
denC = [1 0.25];
C = tf(numC, denC);
L = P*C;
G = feedback(L, 1);
figure;
margin(G);

%The Phase Margin in this plot is show to be 32.6 degrees, which is >= 30
%degrees.

[y, t] = step(G);
steady_state_error_step = abs(1 - y(end));
fprintf('Steady State Error: %f\n', steady_state_error_step);

H1 = freqresp(L, 1);        
tracking_error = abs(1/(1 + H1));
fprintf('Reference Tracking Error: %f\n', tracking_error);

%% 2c)
num = 1;
den = [2 1 0];
P = tf(num, den);
figure;
nyquist(P);
figure;
margin(P);

%These plots are consistent with my sketches in part b.

%% 2d)
G = feedback(P, 1);
figure;
step(G);

%% 2e)

%tau = 0.1
tau = 0.1;
L = tf(num, den, 'InputDelay', tau);
figure;
pzmap(L);
title('Pole-Zero Map (\tau=0.1)');
figure;
nyquist(L);
title('Nyquist Plot (\tau=0.1)');
G = feedback(L, 1);
figure;
step(G);
title('Step Response (\tau=0.1)');

%no. of ccw encirclements around s=-1 = 1
%no. of unstable open-loop poles = 1

%By the Nyquist stability criterion, the closed-loop system should be
%stable, which is what is shown in the step response.

%tau = 0.5
tau = 0.5;
L = tf(num, den, 'InputDelay', tau);
figure;
pzmap(L);
title('Pole-Zero Map (\tau=0.5)');
figure;
nyquist(L);
title('Nyquist Plot (\tau=0.5)');
G = feedback(L, 1);
figure;
step(G);
title('Step Response (\tau=0.5)');

%no. of ccw encirclements around s=-1 = 1
%no. of unstable open-loop poles = 1

%By the Nyquist stability criterion, the closed-loop system should be
%stable, which is what is shown in the step response.

%tau = 1
tau = 1;
L = tf(num, den, 'InputDelay', tau);
figure;
pzmap(L);
title('Pole-Zero Map (\tau=1)');
figure;
nyquist(L);
title('Nyquist Plot (\tau=1)');
G = feedback(L, 1);
figure;
step(G);
title('Step Response (\tau=1)');

%no. of ccw encirclements around s=-1 = 1
%no. of unstable open-loop poles = 1

%By the Nyquist stability criterion, the closed-loop system should be
%stable, which is what is shown in the step response as it approaches the
%steady state.


%tau = 1.5
tau = 1.5;
L = tf(num, den, 'InputDelay', tau);
figure;
pzmap(L);
title('Pole-Zero Map (\tau=1.5)');
figure;
nyquist(L);
title('Nyquist Plot (\tau=1.5)');
G = feedback(L, 1);
figure;
step(G);
title('Step Response (\tau=1.5)');

%no. of ccw encirclements around s=-1 = 0
%no. of unstable open-loop poles = 1

%By the Nyquist stability criterion, the closed-loop system will not be
%stable, which is shown in the step response by its inability to reach
%the steady state value.

%tau = 2
tau = 2;
L = tf(num, den, 'InputDelay', tau);
figure;
pzmap(L);
title('Pole-Zero Map (\tau=2)');
figure;
nyquist(L);
title('Nyquist Plot (\tau=2)');
G = feedback(L, 1);
figure;
step(G);
title('Step Response (\tau=2)');

%no. of ccw encirclements around s=-1 = 0
%no. of unstable open-loop poles = 1

%By the Nyquist stability criterion, the closed-loop system will not be
%stable, which is shown in the step response by its inability to reach
%the steady state value.

%% 3c)
num = 1;
den = [2 1 0];
P = tf(num, den);
P = feedback(1, P);
[mag,~,w] = bode(P);
mag = squeeze(mag);
[peak_mag, idx] = max(mag);
fprintf('Peak Magnitude: %f\n', peak_mag);
figure;
bode(P);

%% 3d)
P_del = P + tf(0.5, [1 1], 'InputDelay', 1);
G_del = feedback(P_del, 1);
figure;
step(G_del);

%The step response validates that the proportional controller can't
%stabilize the true plant at smaller gains there is some underdamping that
%occurs at an earlier time corresponding to the earlier gains.

%% 3e)
figure;
nyquist(P_del);
figure;
pzmap(P_del);

%no. of ccw encirclements around s=-1 = 0
%no. of unstable open-loop poles = 0

%By the Nyquist stability criterion, the closed-loop system will be
%stable. This result is consistent with the step response in part d.