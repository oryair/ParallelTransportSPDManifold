function PlotPositiveMatrix(K)

N      = 20;
vX     = linspace(0, K, N);
vZ     = linspace(0, K, N);
[X, Z] = meshgrid(vX, vZ);
Y      = sqrt(X .* Z);

X = [X, X];
Z = [Z, Z];
Y = [Y, -Y];

% K = convhull(X(:), Y(:), Z(:));
% trisurf(K, X, Y, Z, 'FaceColor', 'r', 'FaceAlpha', 0.05, 'LineStyle', 'None');
plot3(X(:), Y(:), Z(:), '.k');
xlabel('X'); ylabel('Y'); zlabel('Z'); 


end