function f=LE_Chen(t,x)
% LE_Chen.m - Extended system for the Chen chaotic system.
% This function provides the right-hand side for the extended system
% (original system + variational equations) required by FO_Lyapunov.m
% to calculate the Lyapunov exponents.
%
% The Chen system is defined by:
%   dx/dt = a*(y - x)
%   dy/dt = (c - a)*x - x*z + c*y
%   dz/dt = x*y - b*z
%
% Standard chaotic parameters are a=35, b=3, c=28.

% Output data must be a column vector
f=zeros(size(x));

% Parameters for the Chen system
a = 35;
b = 3;
c = 28;

% Variables allocated to the variational equations.
% This is a 3x3 matrix representing the evolution of perturbations.
X= [x(4), x(7), x(10);
    x(5), x(8), x(11);
    x(6), x(9), x(12)];

% Chen system equations
f(1) = a*(x(2) - x(1));
f(2) = (c - a)*x(1) - x(1)*x(3) + c*x(2);
f(3) = x(1)*x(2) - b*x(3);

% Jacobian matrix of the Chen system
J =[-a, a, 0;
    c - a - x(3), c, -x(1);
    x(2), x(1), -b];

% Right-hand side of variational equations (J*X)
f(4:12) = J*X; % To be modified if ne > 3