% --- MATLAB Script to Solve and Plot Fractional-Order Chaotic Systems (Corrected) ---

clear; clc; close all;

% --- 1. System Selection ---
%system_choice = 'Lorenz';
%system_choice = 'Chen';
system_choice = 'Chua';
%system_choice = 'Bhalekar-Gejji';

fprintf('Selected system: %s\n', system_choice);

% --- 2. Define System Parameters and Equations ---
switch system_choice
    case 'Lorenz'
        p.sigma = 10; p.rho = 28; p.beta = 8/3;
        y0 = [0.1, 0.1, 0.1];
        
        % CORRECTED: Use commas for a row vector output
        f = @(t, y, p) [ p.sigma * (y(2) - y(1)), ...
                         y(1) * (p.rho - y(3)) - y(2), ...
                         y(1) * y(2) - p.beta * y(3) ];
        
        plot_title = 'Fractional-Order Lorenz System';

    case 'Chen'
        p.a = 35; p.b = 3; p.c = 28;
        y0 = [-1, 0.5, 1];

        % CORRECTED: Use commas for a row vector output
        f = @(t, y, p) [ p.a * (y(2) - y(1)), ...
                         (p.c - p.a) * y(1) - y(1) * y(3) + p.c * y(2), ...
                         y(1) * y(2) - p.b * y(3) ];
        
        plot_title = 'Fractional-Order Chen System';
    
    %Remark: Does not work for chua, incorrect plot
    case 'Chua'
        p.alpha = 15.6; p.beta = 28; p.m0 = -8/7; p.m1 = -5/7;
        h = @(x, m0, m1) m1*x + 0.5*(m0 - m1)*(abs(x+1) - abs(x-1));
        y0 = [0.1, 0 , 0];

        % CORRECTED: Use commas for a row vector output
        f = @(t, y, p) [ p.alpha * (y(2) - y(1) - h(y(1), p.m0, p.m1)), ...
                         y(1) - y(2) + y(3), ...
                        -p.beta * y(2) ];
        
        plot_title = 'Fractional-Order Chua System';

    case 'Bhalekar-Gejji'
        p.omega = -2.667; p.miu = 10; p.a = 27.3; p.b = 1;
        y0 = [0.7, 0.7, 0.7];

        % CORRECTED: Use commas for a row vector output
        f = @(t, y, p) [ p.omega*y(1)-y(2).^2, ...
                         p.miu*y(3)-p.miu*y(2), ...
                        p.a*y(2)-p.b*y(3)+y(1)*y(2) ];
        
        plot_title = 'Fractional-Order Bhalekar-Geiji System';

    otherwise
        error('Invalid system choice. Check the "system_choice" variable.');
end

% --- 3. Set Solver Parameters ---
alpha = [0.9 0.9 0.9];
tspan = [0, 50];
N = 5000;

% --- 4. Solve the System ---
fprintf('Solving... This may take a moment for high N.\n');
tic;
f_handle = @(t, y) f(t, y, p);
[t, y] = fde_solver_pece_vector(alpha, f_handle, tspan, y0, N);
elapsed_time = toc;
fprintf('Solver finished in %.2f seconds.\n', elapsed_time);

% --- 5. Plot the Results ---
figure('Name', 'Chaotic Attractor', 'NumberTitle', 'off');
plot3(y(:,1), y(:,2), y(:,3), 'b-', 'LineWidth', 0.5);
title(plot_title, 'FontSize', 14);
xlabel('x(t)', 'FontSize', 12);
ylabel('y(t)', 'FontSize', 12);
zlabel('z(t)', 'FontSize', 12);
grid on;
axis tight;
view(-135, 20);