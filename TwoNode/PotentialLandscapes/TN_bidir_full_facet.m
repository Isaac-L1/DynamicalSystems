clear
tic

Xv = linspace(-0.6, 1.4, 301);
Yv = linspace(-0.6, 1.4, 301);
V = zeros(length(Xv), length(Yv));
nus1 = [0.001, 0.005, 0.01, 0.03];
nus2 = [0.001, 0.005, 0.01, 0.03];
beta = 0;
coef = 0.005;
nt = 1000000;
pl = tiledlayout(length(nus1), length(nus2), 'TileSpacing', 'none');
pln = 1;
npe = 3000;

for b = 1:length(nus2)
    
    nu2 = nus2(length(nus2)+1-b);
    
    for a = 1:length(nus1)
        
        nu1 = nus1(a);
        nexttile(pln);
        for i = 1:length(Xv)
            for j = 1:length(Yv)
                y = Yv(j);
                x = Xv(i);
                V(j, i) = (1/4)*x^4 + (1/4)*y^4 - (1/3)*x^3 - (1/3)*y^3 ...
                    - (nu1/2)*x^2 - (nu2/2)*y^2 + nu1*x + nu2*y ...
                    + (beta/2)*x^2 + (beta/2)*y^2 - beta*x*y;
            end
        end
        
        hold on
        contour(Xv, Yv, V, 50)
        %surf(X, Y, V, 'edgeColor', 'none')
        colormap bone
        %xlabel(num2str(nu1))
        %ylabel(num2str(nu2))
        %zlabel('V')
        %title(['\nu 1 = ' num2str(nu1) ', \nu 2 = ' num2str(nu2)])
        %view(3)


        T = 200;
        h = 0.001;
        t = (0:h:T);
        thresh = 0.7;
        th = zeros(1, length(Xv)) + thresh;
        xlb = zeros(1, length(Xv)) + Yv(1);
        xub = zeros(1, length(Xv)) + Yv(length(Yv));
        ylb = zeros(1, length(Xv)) + Xv(1);
        yub = zeros(1, length(Xv)) + Xv(length(Xv));
        kmax = 1;
        alpha = 0.1;
        np = 150;

        plot(Xv, th, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
        plot(th, Yv, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
        plot(Xv, ylb, '-k');
        plot(xlb, Yv, '-k');
        plot(Xv, yub, '-k');
        plot(xub, Yv, '-k');


        for j = 1:kmax %Find kmax realisations with escape times.
            
            rng(1);
            
            nums = randn(4, 4*T/h);
            x = zeros(1, length(t)); %Prepare vector to store realisations
            y = zeros(1, length(t));
            xp = zeros(1, floor(length(t)/np));
            yp = zeros(1, floor(length(t)/np));

            x(1) = -0.1; %Initial conditions
            y(1) = -0.1;
            tx = 0;
            n = 1;
            p=1;
            xesc = 0;
            yesc = 0;
            ne = 0;

            while 1 %Start taking time steps.
                
                if n>=length(nums)
                    nums = randn(4, 4*T/h);
                end
                
                %Evaluate deterministic part using Heun's scheme.
                kx1 = -(x(n)^2 - nu1)*(x(n) - 1) + beta*(y(n)-x(n));
                ky1 = -(y(n)^2 - nu2)*(y(n) - 1) + beta*(x(n)-y(n));

                %Apply noise
                hkx1 = x(n) + h*kx1 + alpha*nums(1, n)*sqrt(h);
                hky1 = y(n) + h*ky1 + alpha*nums(2, n)*sqrt(h);

                kx2 = -(hkx1^2 - nu1)*(hkx1 - 1) + beta*(hky1-hkx1);
                ky2 = -(hky1^2 - nu2)*(hky1 - 1) + beta*(hkx1-hky1);

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
                
                %Check if x has escaped.
                if x(n) > thresh && xesc == 0
                    xesc = 1;
                end

                %Check if y has escaped.
                if y(n) > thresh && yesc == 0
                    yesc = 1;
                end

                %Check if both have escaped.
                if xesc == 1 && yesc == 1
                    ne = ne + 1;
                end
                
                if ne >= npe
                    break;
                end

            end
            %plot(t, x)
            %plot(t, y)
            %figure;
            p = scatter(xp,yp,2,[0,0.4,0.4]);
            %p.Color = [0, 0.6, 0, 1];
        end
        
        
        [roots, vects] = nrfunc(beta, nu1, nu2);
        for i = 1:length(roots)
            switch roots(3, i)
                case 1
                    scatter(roots(1, i), roots(2, i), 'o', 'k');
                case 2
                    scatter(roots(1, i), roots(2, i), 'x', 'k');
                case 3
                    scatter(roots(1, i), roots(2, i), 'd', 'k');
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
        
        
        
        if any([a, (length(nus2)+1-b)] == 1)
            if a == 1
                ly = ylabel(num2str(nus2(length(nus2)+1-b)));
                %ly.Position(1) = ly.Position(1) - 0.05;
            end
            if (length(nus2)+1-b) == 1
                xlabel(num2str(nus1(a)));
            end
        end
        xticklabels({});
        yticklabels({});
        hold off
        pln = pln + 1;
        
    end
end

title(pl, 'Potential landscapes for changing \nu 1 and \nu 2 with unstable manifolds');
xlabel(pl, '\nu 1');
lyt = ylabel(pl, '\nu 2');
lyt.HorizontalAlignment = 'left';

toc
