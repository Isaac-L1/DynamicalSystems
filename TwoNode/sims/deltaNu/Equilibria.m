clear
tic

meanNu = 0.006;
dNs = [-0.005, -0.003, -0.001, 0.001, 0.003, 0.005];
pl = tiledlayout(2, 3, 'TileSpacing', 'compact');
pln = 1;

%One iteration of the loop for each tile.
for deltaNu = 0.001
    nu1 = meanNu-deltaNu;
    nu2 = meanNu+deltaNu;
    sq1 = sqrt(nu1);
    sq2 = sqrt(nu2);
    beta1 = 0.1;
    beta2 = 0;

    nexttile(pln); hold on;
    
    %Parameters for the realisations
    T = 75;
    h = 1e-3;
    t = (0:h:T);
    kmax = 1000;
    alpha = 0.03;
    np = 350;

    
    for j = 1:kmax %Find kmax realisations with escape times.
        nums = randn(4, T/h);
        x = zeros(1, length(t)); %Prepare vector to store realisations
        y = zeros(1, length(t));
        xp = zeros(1, floor(length(t)/np));
        yp = zeros(1, floor(length(t)/np));

        x(1) = -0.1; %Initial conditions
        y(1) = -0.1;
        tx = 0;
        n = 1;
        p=1;

        while n < length(t) %Start taking time steps.
            %Evaluate deterministic part using Heun's scheme.
            kx1 = -(x(n)^2 - nu1)*(x(n) - 1) + beta1*(y(n)-x(n));
            ky1 = -(y(n)^2 - nu2)*(y(n) - 1) + beta2*(x(n)-y(n));

            %Apply noise
            hkx1 = x(n) + h*kx1 + alpha*nums(1, n)*sqrt(h);
            hky1 = y(n) + h*ky1 + alpha*nums(2, n)*sqrt(h);

            kx2 = -(hkx1^2 - nu1)*(hkx1 - 1) + beta1*(hky1-hkx1);
            ky2 = -(hky1^2 - nu2)*(hky1 - 1) + beta2*(hkx1-hky1);

            %Apply heun scheme and noise.
            x(n+1) = x(n) + (h/2)*(kx1+kx2) + alpha*nums(3, n)*sqrt(h);
            y(n+1) = y(n) + (h/2)*(ky1+ky2) + alpha*nums(4, n)*sqrt(h);

            %Step forward
            n = n + 1;
            tx = tx + h;

            if mod(n, np) == 0
                xp(p) = x(n);
                yp(p) = y(n);
                p = p + 1;
            end

        end
        r = plot(xp,yp);
        r.Color = [0, 0.6, 0, 0.2]; %Plot with a low alpha in order to see most likely paths.
    end
    
    roots = nrfuncuni(beta1, beta2, nu1, nu2); %Find the equilibrium points using the 2D Newton-Raphson method.
    for i = 1:length(roots) %Plot each equilibrium with a shape showing stability.
        for j = 1:length(roots)
            switch roots(3, i, j)
                case 1
                    scatter(roots(1, i, j), roots(2, i, j), 'o', 'k');
                case 2
                    scatter(roots(1, i, j), roots(2, i, j), 'x', 'k');
                case 3
                    scatter(roots(1, i, j), roots(2, i, j), 'd', 'k');
            end
        end
    end
    
%     coef = 0.005;
%     nt = 1000000;
%     
%     [roots, vects] = nrfuncuni(beta1, beta2, nu1, nu2);
%     for i = 1:length(roots)
%         switch roots(3, i)
%             case 1
%                 scatter(roots(1, i), roots(2, i), 'o', 'k');
%             case 2
%                 scatter(roots(1, i), roots(2, i), 'x', 'k');
%             case 3
%                 scatter(roots(1, i), roots(2, i), 'd', 'k');
%         end
% 
%         for j = 1:2
%             if vects(3, j, i) == 2
%                 x1 = zeros(1, length(nt)); %Prepare vector to store realisations
%                 y1 = zeros(1, length(nt));
%                 x2 = zeros(1, length(nt));
%                 y2 = zeros(1, length(nt));
% 
%                 x1(1) = roots(1, i) + coef*vects(1, j, i); %Initial conditions
%                 y1(1) = roots(2, i) + coef*vects(2, j, i);
%                 x2(1) = roots(1, i) - coef*vects(1, j, i); %Initial conditions
%                 y2(1) = roots(2, i) - coef*vects(2, j, i);
%                 tx = 0;
%                 n = 1;
% 
%                 while n < nt %Start taking time steps.
% 
%                     kx1 = -(x1(n)^2 - nu1)*(x1(n) - 1) + beta1*(y1(n)-x1(n));
%                     ky1 = -(y1(n)^2 - nu2)*(y1(n) - 1) + beta2*(x1(n)-y1(n));
% 
%                     hkx1 = x1(n) + h*kx1;
%                     hky1 = y1(n) + h*ky1;
% 
%                     kx2 = -(hkx1^2 - nu1)*(hkx1 - 1) + beta1*(hky1-hkx1);
%                     ky2 = -(hky1^2 - nu2)*(hky1 - 1) + beta2*(hkx1-hky1);
%                     
%                     x1(n+1) = x1(n) + (h/2)*(kx1+kx2);
%                     y1(n+1) = y1(n) + (h/2)*(ky1+ky2);
% 
%                     %Step forward
%                     n = n + 1;
% 
%                 end
% 
%                 n = 1;
% 
%                 while n < nt %Start taking time steps.
% 
%                     %Evaluate deterministic part using Heun's scheme.
%                     kx1 = -(x2(n)^2 - nu1)*(x2(n) - 1) + beta1*(y2(n)-x2(n));
%                     ky1 = -(y2(n)^2 - nu2)*(y2(n) - 1) + beta2*(x2(n)-y2(n));
% 
%                     hkx1 = x2(n) + h*kx1;
%                     hky1 = y2(n) + h*ky1;
% 
%                     kx2 = -(hkx1^2 - nu1)*(hkx1 - 1) + beta1*(hky1-hkx1);
%                     ky2 = -(hky1^2 - nu2)*(hky1 - 1) + beta2*(hkx1-hky1);
% 
%                     %Apply heun scheme and noise.
%                     x2(n+1) = x2(n) + (h/2)*(kx1+kx2);
%                     y2(n+1) = y2(n) + (h/2)*(ky1+ky2);
% 
%                     %Step forward
%                     n = n + 1;
% 
%                 end
%                 %plot(t, x)
%                 %plot(t, y)
%                 %figure;
%                 plot(x1,y1,'r');
%                 plot(x2,y2,'r');
%             end
%         end
%     end
    
    title(['\delta \nu = ' num2str(deltaNu)]);
    xlabel('X1');
    ylabel('X2');
    ax = linspace(-0.5, 1.5, 10);
    r = zeros(1,10);
    
%     plot(ax,r-sq2);
%     plot(ax,r+sq2);
%     plot(ax,r+1);
%     plot(r-sq1,ax);
%     plot(r+sq1,ax);
%     plot(r+1,ax);
    
    hold off
    pln = pln + 1;
end

title(pl, 'Mean \nu = 0.006 realisations and equilibria');
toc



