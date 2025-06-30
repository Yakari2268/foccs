% Fractional order
alpha = 0.9;

% System function: f(t, y) returns a column vector
f = @(t, y) [y(2); -y(1)];

% Time settings
t0 = 0;
T = 10;
h = 0.01;

% Initial conditions (column vector)
y0 = [0; 1];

% Call fde12
[t, y] = fde12(alpha, f, t0, T, y0, h);

% Plot results
figure;
plot(t, y(1,:), 'b', t, y(2,:), 'r--');
legend('y_1(t)', 'y_2(t)');
xlabel('t'); ylabel('y(t)');
title('Fractional System: D^{0.9}y_1 = y_2, D^{0.9}y_2 = -y_1');
grid on;
