clear

%% Variables

alpha = 0.1;
betas = [0,0,0;0,0,0;0,0,0];
nu1 = 0.01;
nu2 = 0.01;
nu3 = 0.01;

T = 100;
h = 1e-3;
t = (0:h:T);

kmax = 1;

%% Simulation
tic

thresh = 0.6;
th = zeros(1, length(t)) + thresh;

figure; hold on

for j = 1:kmax %Find kmax realisations with escape times.
    nums = randn(6, T/h);
    x1 = zeros(1, length(t)); %Prepare vector to store realisations
    x2 = zeros(1, length(t));
    x3 = zeros(1, length(t));
    
    x1(1) = -0.1; %Initial conditions
    x2(1) = -0.1;
    x3(1) = -0.1;
    tx = 0;
    n = 1;
        
    while n < length(t) %Start taking time steps.
        %Evaluate deterministic part using Heun's scheme.
        kx11 = -(x1(n)^2 - nu1)*(x1(n) - 1) + betas(1, 2)*(x2(n)-x1(n)) + betas(1, 3)*(x3(n)-x1(n));
        kx21 = -(x2(n)^2 - nu2)*(x2(n) - 1) + betas(2, 1)*(x1(n)-x2(n)) + betas(2, 3)*(x3(n)-x2(n));
        kx31 = -(x3(n)^2 - nu3)*(x3(n) - 1) + betas(3, 1)*(x1(n)-x3(n)) + betas(3, 2)*(x2(n)-x3(n));

        %Apply noise
        hkx11 = x1(n) + h*kx11 + alpha*nums(1, n)*sqrt(h);
        hkx21 = x2(n) + h*kx21 + alpha*nums(2, n)*sqrt(h);
        hkx31 = x3(n) + h*kx31 + alpha*nums(3, n)*sqrt(h);

        kx12 = -(hkx11^2 - nu1)*(hkx11 - 1) + betas(1, 2)*(hkx21-hkx11) + betas(1, 3)*(hkx31-hkx11);
        kx22 = -(hkx21^2 - nu2)*(hkx21 - 1) + betas(2, 1)*(hkx11-hkx21) + betas(2, 3)*(hkx31-hkx21);
        kx32 = -(hkx31^2 - nu3)*(hkx31 - 1) + betas(3, 1)*(hkx11-hkx31) + betas(3, 2)*(hkx21-hkx31);

        %Apply heun scheme and noise.
        x1(n+1) = x1(n) + (h/2)*(kx11+kx12) + alpha*nums(4, n)*sqrt(h);
        x2(n+1) = x2(n) + (h/2)*(kx21+kx22) + alpha*nums(5, n)*sqrt(h);
        x3(n+1) = x3(n) + (h/2)*(kx31+kx32) + alpha*nums(6, n)*sqrt(h);

        %Step forward
        n = n + 1;
        tx = tx + h;
        
    end
    plot(t, x1)
    plot(t, x2)
    plot(t, x3)
    
end

%plot(t, th)
%hold off

toc
