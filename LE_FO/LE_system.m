

% Parameters for the Bhalekar-Gejji system
p.omega = -2.667;
p.mu = 10;
p.a = 27.3;
p.b = 1;
p.k = 1;
p.k1 = 26.79;

% Bhalekar-Gejji system equations
f = @(t,x) [p.omega*x(1) - x(2)^2 ,...
            p.mu*(x(3) - x(2)) ,...
            p.a*x(2) - p.b*x(3) + x(1)*x(2)*(1-p.k*(sin(p.k1*x(2))))];

ne = 3;
q = 0.98;
n = 5000;
x0 = [0.7 0.7 0.7];
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

q = 0.98;ne = 3;
%--- Find the lyapunov exponent---
fprintf("Finding the lyapunov exponent");

[t_hist,LE_hist] = FO_Lyapunov(ne,@mw_Bhalekar_Gejji,0,0.02,3000, [0.1;0.1;0.1], 0.005,q, 1000);
plot_LE(t_hist,LE_hist,q,ne)
