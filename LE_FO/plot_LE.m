function [] = plot_LE(t,LE,q,ne)
% --- Plot the Results ---
disp('Plotting results...');
figure;
plot(t, LE, 'LineWidth', 1.5);

% --- Add plot decorations ---
title(['Evolution of Lyapunov Exponents (q = ' num2str(q) ')']);
xlabel('Time (t)');
ylabel('Lyapunov Exponents (\lambda_i)');
grid on;
ax = gca;
ax.FontSize = 12;

% Create a legend
legend_labels = cell(ne, 1);
for i = 1:ne
    legend_labels{i} = ['$\lambda_{' num2str(i) '}$'];
end
legend(legend_labels, 'Interpreter', 'latex', 'Location', 'northeast');

disp('Done.');
end