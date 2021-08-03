clear

%% Variables

alpha = 0.1;
%betas = logspace(-4, -0.7, 50);
%nus = logspace(-3, -1.5, 50);
beta  = 0;
nus1 = logspace(-3, -1.5, 80);
nus2 = logspace(-3, -1.5, 80);

kmax = 2000;

T = 1000;
h = 1e-3;

Ex = zeros(length(nus1), length(nus2));   % holder for mean escape times
Ey = zeros(length(nus1), length(nus2));
Emax = zeros(length(nus1), length(nus2));
Emin = zeros(length(nus1), length(nus2));
D = zeros(length(nus1), length(nus2));
R = zeros(length(nus1), length(nus2));

%% Simulation
tic
thresh = 0.8;

for b = 1:length(nus1)
    
    nu1 = nus1(b);
    fprintf(['\n Computing escape times for nu1 = ' num2str(nu1) '\n']); %Print which step it is on

    parfor i = 1:length(nus2)

        nu2 = nus2(i);
        fprintf(['\n Computing escape times for nu2 = ' num2str(nu2) '\n']); %Print which step it is on

        Exn = zeros(1, kmax); %Set up escape time vector
        Eyn = zeros(1, kmax); %Set up escape time vector

        for j = 1:kmax %Find kmax realisations with escape times.

            x = zeros(1, T/h); %Prepare vector to store realisations
            y = zeros(1, T/h);

            x(1) = -0.1; %Initial conditions
            y(1) = -0.1;
            tx = 0;
            n = 1;
            xesc = 0;
            yesc = 0;
            
            while 1 %Start taking time steps.
                
                %Evaluate deterministic part using Heun's scheme.
                kx1 = F(x(n), y(n), nu1, beta);
                ky1 = F(y(n), x(n), nu2, beta);

                %Apply noise
                hkx1 = h*kx1 + alpha*randn*sqrt(h);
                hky1 = h*ky1 + alpha*randn*sqrt(h);

                kx2 = F((x(n) + hkx1), (y(n)+ hky1), nu1, beta);
                ky2 = F((y(n) + hky1), (x(n)+ hkx1), nu2, beta);

                %Apply heun scheme and noise.
                x(n+1) = x(n) + (h/2)*(kx1+kx2) + alpha*randn*sqrt(h);
                y(n+1) = y(n) + (h/2)*(ky1+ky2) + alpha*randn*sqrt(h);
                
                %Step forward
                n = n + 1;
                tx = tx + h;

                %Check if x has escaped.
                if x(n) > thresh && xesc == 0
                    Exn(j) = tx;
                    xesc = 1;
                end

                %Check if y has escaped.
                if y(n) > thresh && yesc == 0
                    Eyn(j) = tx;
                    yesc = 1;
                end

                %Check if both have escaped.
                if xesc == 1 && yesc == 1
                    break
                end

                
            end
        end

        Ex(b, i) = mean(Exn);
        Ey(b, i) = mean(Eyn);
        Emax(b, i) = max([Ex(b, i), Ey(b, i)]);
        Emin(b, i) = min([Ex(b, i), Ey(b, i)]);
        D(b, i) = (Emax(b, i) - Emin(b, i));
        R(b, i) = (D(b, i) / Emin(b, i));

    end
end
toc

%% Plot the means 3d

figure;hold on;
surf(nus1, nus2, R, 'edgeColor', 'none')
colormap turbo
% hleg = legend('0', '0.002', '0.02', '0.2', 'Location', 'NW');
% htitle = get(hleg,'Title');
% set(htitle,'String','Beta')
xlabel('\nu 1')
ylabel('\nu 2','Rotation',0)
zlabel('E(\tau)')
title('Mean escape times for two nodes coupled bidirectionally against \nu')
set(gca,'xscale','log');
set(gca,'yscale','log');
view(3)
hold off

%% Plot the means 2d

figure;hold on;
plot(nus1, Ex);
colormap turbo
% hleg = legend('0', '0.002', '0.02', '0.2', 'Location', 'NW');
% htitle = get(hleg,'Title');
% set(htitle,'String','Beta')
xlabel('\nu')
ylabel('E(\tau)','Rotation',0)
title('Mean escape times for two nodes coupled bidirectionally against \nu')
set(gca,'xscale','log');
hold off

%% Plot the sds
figure;hold on;
plot(nus, Sx);
plot(nus, Sy);

hleg = legend('x', 'y', 'Location', 'NW');
htitle = get(hleg,'Title');
set(htitle,'String','Node')
xlabel('\nu')
ylabel('E[\tau]','Rotation',0)
title('SD of escape times for two nodes coupled bidirectionally against \nu')
set(gca,'xscale','log');
set(gca,'yscale','log');
box on
hold off



%% Function Declarations

function [x1] = F(x0, y0, v, B)
    
    x1 = -(x0^2 - v)*(x0 - 1) + B*(y0-x0);

end
