function f = mw_ext_fun(t,x)
% Output data must be a column vector
f=zeros(size(x));

% Parameters for the multi-wing chaotic system
a = 35;
b = 3;
c = 28;
k = 0.5;
k1 = 1;

% Variables allocated to the variational equations (perturbation matrix)
X= [x(4), x(7), x(10);
    x(5), x(8), x(11);
    x(6), x(9), x(12)];

% Bhalekar-Gejji system equations
f(1) = a*(x(2) - x(1));
f(2) = c*x(2) - x(1)*x(3)*(1-k*sin(k1*x(3))) + (c-a)*x(1);
f(3) = -b*x(3) + x(1)*x(2);

% Jacobian matrix of the Bhalekar-Gejji system
J = [-a, a, 0;
     -x(3)*(1-k*sin(k1*x(3))) + (c-a), c, -x(1) + k*x(1)*(sin(k1*x(3)) + k1*x(3)*cos(k1*x(3)));
     x(2), x(1), -b];

% Right-hand side of variational equations (J*X)
f(4:12) = J*X; % To be modified if ne > 3

end