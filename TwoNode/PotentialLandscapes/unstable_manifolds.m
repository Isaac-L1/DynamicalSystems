clear

beta = 0;
nu1 = 0.01;
nu2 = 0.01;
coef = 0.005;
nt = 1000000;
h = 0.001;

X = linspace(-0.6, 1.4, 301);
Y = linspace(-0.6, 1.4, 301);
V = zeros(length(X), length(Y));

figure; hold on;

for i = 1:length(X)
    for j = 1:length(Y)
        y = Y(j);
        x = X(i);
        V(j, i) = (1/4)*x^4 + (1/4)*y^4 - (1/3)*x^3 - (1/3)*y^3 ...
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

[roots, vects] = nrfunc(beta, nu1, nu2);

for i = 1:length(roots)
    switch roots(3, i)
        case 1
            scatter(roots(2, i), roots(1, i), 'o', 'k');
        case 2
            scatter(roots(2, i), roots(1, i), 'x', 'k');
        case 3
            scatter(roots(2, i), roots(1, i), 'd', 'k');
    end
    
    for j = 1:2
        if vects(3, j, i) == 2
            x1 = zeros(1, length(nt)); %Prepare vector to store realisations
            y1 = zeros(1, length(nt));
            x2 = zeros(1, length(nt));
            y2 = zeros(1, length(nt));

            x1(1) = roots(1, i) + coef*vects(1, j, i); %Initial conditions
            y1(1) = roots(2, i) + coef*vects(2, j, i);
            x2(1) = roots(1, i) - coef*vects(1, j, i); %Initial conditions
            y2(1) = roots(2, i) - coef*vects(2, j, i);
            tx = 0;
            n = 1;

            while n < nt %Start taking time steps.
                
                %Evaluate deterministic part using Heun's scheme.
                kx1 = -(x1(n)^2 - nu1)*(x1(n) - 1) + beta*(y1(n)-x1(n));
                ky1 = -(y1(n)^2 - nu2)*(y1(n) - 1) + beta*(x1(n)-y1(n));
                
                hkx1 = x1(n) + h*kx1;
                hky1 = y1(n) + h*ky1;
                
                kx2 = -(hkx1^2 - nu1)*(hkx1 - 1) + beta*(hky1-hkx1);
                ky2 = -(hky1^2 - nu2)*(hky1 - 1) + beta*(hkx1-hky1);

                %Apply heun scheme and noise.
                x1(n+1) = x1(n) + (h/2)*(kx1+kx2);
                y1(n+1) = y1(n) + (h/2)*(ky1+ky2);
                
                %Step forward
                n = n + 1;
                
            end
            
            n = 1;
            
            while n < nt %Start taking time steps.
                
                %Evaluate deterministic part using Heun's scheme.
                kx1 = -(x2(n)^2 - nu1)*(x2(n) - 1) + beta*(y2(n)-x2(n));
                ky1 = -(y2(n)^2 - nu2)*(y2(n) - 1) + beta*(x2(n)-y2(n));
                
                hkx1 = x2(n) + h*kx1;
                hky1 = y2(n) + h*ky1;
                
                kx2 = -(hkx1^2 - nu1)*(hkx1 - 1) + beta*(hky1-hkx1);
                ky2 = -(hky1^2 - nu2)*(hky1 - 1) + beta*(hkx1-hky1);

                %Apply heun scheme and noise.
                x2(n+1) = x2(n) + (h/2)*(kx1+kx2);
                y2(n+1) = y2(n) + (h/2)*(ky1+ky2);
                
                %Step forward
                n = n + 1;
                
            end
            %plot(t, x)
            %plot(t, y)
            %figure;
            plot(x1,y1,'r');
            plot(x2,y2,'r');
        end
    end
end

