[t,x] = fde12(0.9,@chen_system,0,300,[1;1;1],0.02);

figure;
plot3(x(:,1), x(:,2), x(:,3));
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
title('Chen Attractor (Fractional Order)');
grid on;