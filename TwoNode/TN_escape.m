clear

%% Variables

alpha = 0.1;
% betas = logspace(-4, -0.7, 50);
nus1 = logspace(-3, -1.5, 10);
nus2 = logspace(-3, -1.5, 10);
%nus1 = logspace(-3, -1.5, 5);
%nus2 = logspace(-3, -1.5, 5);
beta = 0;
% nus1 = 0.01;
% nus2 = 0.01;

kmax = 10;

T = 250;
h = 1e-3;

Ex = zeros(length(nus1), length(nus2));   % holder for mean escape times
Ey = zeros(length(nus1), length(nus2));
Sx = zeros(length(nus1), length(nus2));   % holder for stdevs
Sy = zeros(length(nus1), length(nus2));

%% Simulation
thresh = 0.8;
tic

for b = 1:length(nus2)
    
    nu2 = nus2(b);
    %fprintf(['\n Computing escape times for beta = ' num2str(beta) '\n']); %Print which step it is on

    for i = 1:length(nus1)

        nu1 = nus1(i);
        %fprintf(['\n Computing escape times for nu = ' num2str(nu1) '\n']); %Print which step it is on

        Exn = zeros(1, kmax); %Set up escape time vector
        Eyn = zeros(1, kmax); %Set up escape time vector
%         needRegen = 0;
%         if kmax <= 300
%             nums = randn(4, kmax*T/h);
%         else
%             nums = randn(4, 300*T/h);
%             needRegen = 1;
%         end

        for j = 1:kmax %Find kmax realisations with escape times.
            
            nums = randn(4, 2*T/h);
            
            x = zeros(1, T/h); %Prepare vector to store realisations
            y = zeros(1, T/h);

            x(1) = -0.1; %Initial conditions
            y(1) = -0.1;
            tx = 0;
            n = 1;
            xesc = 0;
            yesc = 0;
            
            
            while 1 %Start taking time steps.
%                 nn = (mod(j+299, 300))*(T/h)+n;
                
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
            
%             if needRegen && mod(j, 300) == 299
%                 nums = randn(4, 300*T/h);
%             end
        end

        Ex(b, i) = mean(Exn);
        Ey(b, i) = mean(Eyn);
        Sx(b, i) = std(Exn);
        Sy(b, i) = std(Eyn);

    end
end
toc

% %% Plot the means 3d
% 
% figure;hold on;
% surf(betas, nus, Ex, 'edgeColor', 'none')
% colormap turbo
% % hleg = legend('0', '0.002', '0.02', '0.2', 'Location', 'NW');
% % htitle = get(hleg,'Title');
% % set(htitle,'String','Beta')
% xlabel('\beta')
% ylabel('\nu','Rotation',0)
% zlabel('E(\tau)')
% title('Mean escape times for two nodes coupled bidirectionally against \nu')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% view(3)
% hold off
% 
% %% Plot the means 2d
% 
% figure;hold on;
% plot(nus1, E);
% colormap turbo
% % hleg = legend('0', '0.002', '0.02', '0.2', 'Location', 'NW');
% % htitle = get(hleg,'Title');
% % set(htitle,'String','Beta')
% xlabel('\nu')
% ylabel('E(\tau)','Rotation',0)
% title('Mean escape times for two nodes coupled bidirectionally against \nu')
% set(gca,'xscale','log');
% hold off
% 
% %% Plot the sds
% figure;hold on;
% plot(nus, Sx);
% plot(nus, Sy);
% 
% hleg = legend('x', 'y', 'Location', 'NW');
% htitle = get(hleg,'Title');
% set(htitle,'String','Node')
% xlabel('\nu')
% ylabel('E[\tau]','Rotation',0)
% title('SD of escape times for two nodes coupled bidirectionally against \nu')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% box on
% hold off
% 
% 
% 
%% Function Declarations
% 
% function [x1] = F(x0, y0, v, B)
%     
%     x1 = -(x0^2 - v)*(x0 - 1) + B*(y0-x0);
% 
% end

