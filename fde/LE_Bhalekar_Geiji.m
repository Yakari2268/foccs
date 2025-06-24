function f = LE_Bhalekar_Geiji(t,x)
 p.omega = -2.667; p.miu = 10; p.a = 27.3; p.b = 1;

 X= [x(4), x(7), x(10);
    x(5), x(8), x(11);
    x(6), x(9), x(12)];
 % CORRECTED: Use commas for a row vector output
 f = @(x, p) [ p.omega*x(1)-x(2).^2;
                  p.miu*x(3)-p.miu*x(2);
                 p.a*x(2)-p.b*x(3)+x(1)*x(2)];

 J = [p.omega, -2*x(2), 0;
     0, -p.miu, p.miu;
     x(2), p.a+x(1), -p.b];

 f(4:12) = J*X; % To be modified if ne > 3
 