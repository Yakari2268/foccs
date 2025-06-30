
clear; clc; close all;
%Parameters and their value
a = 35; b = 3; c = 28; k = 0.5; k1 = 1;


f = @(t,x) [a*(x(2) - x(1)), ...
            c*x(2) - x(1)*x(3)*(1-k*sin(k1*(x(3)))) + (c-a)*x(1),...
            -b*x(3) + x(1)*x(2)];
ne = 3;
q = 0.98;
n = 0.0004;
x0 = [-1 0.5 1];
f_handle = @(t,x) f(t,x);

% --- 4. Solve the System ---
fprintf('Solving... This may take a moment for high N.\n');
tic
[t,y] = fde_solver_pece_vector(q,f_handle,[0 50],x0,n);
t_elapsed = toc;

fprintf("Time elapsed = %.2f", t_elapsed);

figure(1)
plot3(y(:,1),y(:,2),y(:,3),'b-','LineWidth',0.5);
grid on
axis tight

%--- Find the lyapunov exponent---
fprintf("Finding the lyapunov exponent");

[t_hist,LE_hist] = FO_Lyapunov(ne,@mw_ext_fun,0,0.02,300, [0.1;0.1;0.1], 0.005,q, 1000);
plot_LE(t_hist,LE_hist,q,ne)