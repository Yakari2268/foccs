function f = LE_Bhalekar_Gejji(t,x)
% LE_Bhalekar_Gejji.m - Extended system for the Bhalekar-Gejji chaotic system.
% This function provides the right-hand side for the extended system
% (original system + variational equations) required by FO_Lyapunov.m.
%
% The Bhalekar-Gejji system is defined by:
%   dx/dt = omega*x - y^2
%   dy/dt = mu*(z - y)
%   dz/dt = a*y - b*z + x*y
%
% Standard chaotic parameters are omega=-2.667, mu=10, a=27.3, b=1.

% Output data must be a column vector
f=zeros(size(x));

% Parameters for the Bhalekar-Gejji system
p.omega = -2.667;
p.mu = 10;
p.a = 27.3;
p.b = 1;

% Variables allocated to the variational equations (perturbation matrix)
X= [x(4), x(7), x(10);
    x(5), x(8), x(11);
    x(6), x(9), x(12)];

% Bhalekar-Gejji system equations
f(1) = p.omega*x(1) - x(2)^2;
f(2) = p.mu*(x(3) - x(2));
f(3) = p.a*x(2) - p.b*x(3) + x(1)*x(2);

% Jacobian matrix of the Bhalekar-Gejji system
J = [p.omega, -2*x(2), 0;
     0, -p.mu, p.mu;
     x(2), p.a+x(1), -p.b];

% Right-hand side of variational equations (J*X)
f(4:12) = J*X; % To be modified if ne > 3