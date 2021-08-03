clear
tic
X = linspace(-0.6, 1.4, 301);
Y = linspace(-0.6, 1.4, 301);
V = zeros(length(X), length(Y));
nu1 = 0.00291962;
nu2 = 0.00291962;
beta = 0.00644236;
figure; hold on;

for i = 1:length(X)
    for j = 1:length(Y)
        y = Y(j);
        x = X(i);
        V(i, j) = (1/4)*x^4 + (1/4)*y^4 - (1/3)*x^3 - (1/3)*y^3 ...
            - (nu1/2)*x^2 - (nu2/2)*y^2 + nu1*x + nu2*y ...
            + (beta/2)*x^2 + (beta/2)*y^2 - beta*x*y;
    end
end

contour(X, Y, V, 50)
%surf(X, Y, V, 'edgeColor', 'none')
colormap turbo
xlabel('X')
ylabel('Y')
%zlabel('V')
title('Surface showing the Potential Landscape')
%view(3)

roots = nrfunc(beta, nu1, nu2);
for i = 1:length(roots)
    switch roots(3, i)
        case 1
            scatter(roots(2, i), roots(1, i), 'o', 'k');
        case 2
            scatter(roots(2, i), roots(1, i), 'x', 'k');
        case 3
            scatter(roots(2, i), roots(1, i), 'd', 'k');
    end
end


T = 75;
h = 1e-3;
t = (0:h:T);
thresh = 0.7;
th = zeros(1, length(X)) + thresh;
kmax = 1002;
alpha = 0.1;
np = 350;

plot(X, th, '-k');
plot(th, Y, '-k');


% for j = 1:kmax %Find kmax realisations with escape times.
%     nums = randn(4, T/h);
%     x = zeros(1, length(t)); %Prepare vector to store realisations
%     y = zeros(1, length(t));
%     xp = zeros(1, floor(length(t)/np));
%     yp = zeros(1, floor(length(t)/np));
% 
%     x(1) = -0.1; %Initial conditions
%     y(1) = -0.1;
%     tx = 0;
%     n = 1;
%     p=1;
% 
%     while n < length(t) %Start taking time steps.
%         %Evaluate deterministic part using Heun's scheme.
%         kx1 = -(x(n)^2 - nu1)*(x(n) - 1) + beta*(y(n)-x(n));
%         ky1 = -(y(n)^2 - nu2)*(y(n) - 1) + beta*(x(n)-y(n));
% 
%         %Apply noise
%         hkx1 = x(n) + h*kx1 + alpha*nums(1, n)*sqrt(h);
%         hky1 = y(n) + h*ky1 + alpha*nums(2, n)*sqrt(h);
% 
%         kx2 = -(hkx1^2 - nu1)*(hkx1 - 1) + beta*(hky1-hkx1);
%         ky2 = -(hky1^2 - nu2)*(hky1 - 1) + beta*(hkx1-hky1);
% 
%         %Apply heun scheme and noise.
%         x(n+1) = x(n) + (h/2)*(kx1+kx2) + alpha*nums(3, n)*sqrt(h);
%         y(n+1) = y(n) + (h/2)*(ky1+ky2) + alpha*nums(4, n)*sqrt(h);
% 
%         %Step forward
%         n = n + 1;
%         tx = tx + h;
% 
%         if mod(n, np) == 0
%             xp(p) = x(n);
%             yp(p) = y(n);
%             p = p + 1;
%         end
% 
%     end
%     %plot(t, x)
%     %plot(t, y)
%     %figure;
%     r = plot(xp,yp);
%     r.Color = [0, 0.6, 0, 0.1];
% end
toc