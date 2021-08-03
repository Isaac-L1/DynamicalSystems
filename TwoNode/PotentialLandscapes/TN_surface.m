clear

%% Set up
X = linspace(-0.6, 1.4, 301);
Y = linspace(-0.6, 1.4, 301);
V = zeros(length(X), length(Y));
nu1 = 0.01;
nu2 = 0.01;
beta = 0;

%% Find potential

for i = 1:length(X)
    for j = 1:length(Y)
        y = Y(j);
        x = X(i);
        V(i, j) = (1/4)*x^4 + (1/4)*y^4 - (1/3)*x^3 - (1/3)*y^3 ...
            - (nu1/2)*x^2 - (nu2/2)*y^2 + nu1*x + nu2*y ...
            + (beta/2)*x^2 + (beta/2)*y^2 - beta*x*y;
    end
end

%% Plot it
figure; hold on;
surf(X, Y, V, 'edgeColor', 'none')
colormap turbo
xlabel('X')
ylabel('Y')
zlabel('V')
title('Surface showing the Potential Landscape')
view(3)