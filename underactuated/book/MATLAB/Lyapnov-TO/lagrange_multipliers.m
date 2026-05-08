% cost function
J = @(x1, x2) 2*x1.^2 - 8*x1 - 4*x2 + 4*x2.^2 + 12;

% Optimal point
x1_opt = 5/3;
x2_opt = 1/3;

[x1, x2] = meshgrid(-2:0.01:4, -2:0.01:4);
Z = J(x1, x2);
J(x1_opt, x2_opt)

figure;
contour(x1, x2, Z, 80, 'LineWidth', 1.2);  % 80 contour levels
hold on;

% constraint line x1 + x2 = 2
x1_line = linspace(-2, 4, 100);
x2_line = 2 - x1_line;
plot(x1_line, x2_line, 'r', 'LineWidth', 2);

% optimal point
plot(x1_opt, x2_opt, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);

xlabel('x_1');
ylabel('x_2');
title('Cost Function Contour Plot with Constraint');
legend('Cost Contours', 'Constraint: x_1 + x_2 = 2', 'Optimal Point', 'Location', 'northwest');
grid on;
axis equal;
colorbar;
hold off;
