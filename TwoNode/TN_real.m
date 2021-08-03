clear

%% Variables

alpha = 0.1;
beta = 0.002;
nu = 0.01;

T = 100;
h = 1e-3;
t = (0:h:T);

kmax = 5;

%% Simulation
tic

thresh = 0.6;
th = zeros(1, length(t));
th(1) = thresh;

figure; hold on

for j = 1:kmax %Find kmax realisations with escape times.
    nums = randn(4, T/h);
    x = zeros(1, length(t)); %Prepare vector to store realisations
    y = zeros(1, length(t));
    
    x(1) = -0.1; %Initial conditions
    y(1) = -0.1;
    tx = 0;
    n = 1;
        
    while n < length(t) %Start taking time steps.
        %Evaluate deterministic part using Heun's scheme.
        kx1 = F(x(n), y(n), nu, beta);
        ky1 = F(y(n), x(n), nu, beta);
        
        hkx1 = h*kx1 + alpha*nums(1, n)*sqrt(h);
        hky1 = h*ky1 + alpha*nums(2, n)*sqrt(h);
        
        kx2 = F((x(n) + hkx1), (y(n)+ hky1), nu, beta);
        ky2 = F((y(n) + hky1), (x(n)+ hkx1), nu, beta);
        
        x(n+1) = x(n) + (h/2)*(kx1+kx2) + alpha*nums(3, n)*sqrt(h);
        y(n+1) = y(n) + (h/2)*(ky1+ky2) + alpha*nums(4, n)*sqrt(h);
            
        th(n+1) = thresh;
            
        tx = tx + h;
            
            
        %Step forward
        n = n + 1;
    end
    %plot(t, x)
    %plot(t, y)
    %figure;
    plot(x,y);
end

%plot(t, th)
%hold off

toc

%% Function Declarations

function [x1] = F(x0, y0, v, B)
    
    x1 = -(x0^2 - v)*(x0 - 1) + B*(y0-x0);

end
