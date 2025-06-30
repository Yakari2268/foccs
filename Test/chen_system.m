function dxdt = chen_system(t,x)
    % Chen System equations
    %   dx/dt = a*(y - x)
    %   dy/dt = (c - a)*x - x*z + c*y
    %   dz/dt = x*y - b*z

    % Parameters
    a = 35; b = 3; c = 28;

    % Derivatives
    dxdt = zeros(3,1);
    dxdt(1) = a * (x(2) - x(1));
    dxdt(2) = (c - a) * x(1) - x(1) * x(3) + c * x(2);
    dxdt(3) = x(1) * x(2) - b * x(3);
end
