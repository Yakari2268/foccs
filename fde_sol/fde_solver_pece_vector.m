function [t, y] = fde_solver_pece_vector(alpha, f, tspan, y0, N)
% Solves a system of FDEs D^alpha y = f(t,y) for 0 < alpha < 1 using the
% Adams-Bashforth-Moulton Predictor-Corrector method.
%
% INPUTS:
%   alpha - Fractional order. Can be a scalar (if all equations have the
%           same order) or a row vector (one order per equation).
%   f     - Function handle f(t,y), where y is a row vector. f must
%           return a row vector of the same dimension as y.
%   tspan - Time interval [t_start, t_end].
%   y0    - Row vector of initial conditions [y1_0, y2_0, ...].
%   N     - Number of time steps to take.
%
% OUTPUTS:
%   t     - Column vector of time points.
%   y     - Solution matrix where each row corresponds to a time point and
%           each column corresponds to a dependent variable.

% --- 1. SETUP ---
t_start = tspan(1);
t_end = tspan(2);
h = (t_end - t_start) / N; % Calculate step size
t = linspace(t_start, t_end, N + 1)'; % Create time grid as a column vector

% Determine the number of equations from the initial condition vector
num_eqns = length(y0);

% If alpha is a scalar, expand it to a vector
if isscalar(alpha)
    alpha = alpha * ones(1, num_eqns);
end

% Initialize solution matrix and set initial conditions
y = zeros(N + 1, num_eqns);
y(1, :) = y0;

% Store history of f(t,y) values for the convolution sums
f_history = zeros(N + 1, num_eqns);
f_history(1, :) = f(t(1), y(1, :));

% Pre-calculate gamma terms (these are vectors now)
for i=1:length(alpha)
    ga1(i) = gamma_lanczos(alpha(i) + 1);
    ga2(i) = gamma_lanczos(alpha(i) + 2);
end

% --- 2. MAIN LOOP (Vectorized) ---
% Iterate from k=1 to N, calculating the solution at each step t(k+1)
for k = 1:N
    
    % --- PREDICTOR STEP ---
    predictor_sum = zeros(1, num_eqns);
    for j = 0:k-1
        % Adams-Bashforth weights (b_j) are now vectors
        b_j = (k-j).^alpha - (k-1-j).^alpha;
        predictor_sum = predictor_sum + b_j .* f_history(j+1, :);
    end
    
    y_p = y0 + (h.^alpha ./ ga1) .* predictor_sum;
    
    % Evaluate f at the predicted point
    f_p = f(t(k+1), y_p);
    
    % --- CORRECTOR STEP ---
    corrector_sum = zeros(1, num_eqns);
    
    % Contribution from j=0
    a0 = (k-1).^(alpha+1) - (k-1-alpha) .* k.^alpha;
    corrector_sum = corrector_sum + a0 .* f_history(1, :);
    
    % Contribution from j=1 to k-1
    for j = 1:k-1
        % Adams-Moulton weights (a_j) are now vectors
        aj = (k-j+1).^(alpha+1) + (k-j-1).^(alpha+1) - 2*(k-j).^(alpha+1);
        corrector_sum = corrector_sum + aj .* f_history(j+1, :);
    end
    
    % Combine terms using element-wise operations for vectors
    y(k+1, :) = y0 + (h.^alpha ./ ga2) .* (f_p + corrector_sum);
    
    % Update history with the f value from the new corrected solution
    f_history(k+1, :) = f(t(k+1), y(k+1, :));
end

end