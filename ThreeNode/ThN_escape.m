clear

N = 3;
alpha = 0.1;
betas = [0,0,0;0,0,0;0,0,0];
nu1 = 0.01;
nu2 = 0.01;
nu3 = 0.01;
h = 1e-3;
sh = sqrt(h);
kmax = 1000;
thresh = 0.6;

%%
tic

%Set up escape time matrices and maximum time.
Ex = zeros(N, 1);
Exn = zeros(N, kmax);

%Iterate through number of realisations.
for kn = 1:kmax
    
    nums = randn(6, 500000);
    %Set up vectors.
    tx = 0;
    n = 1;
    x1 = -0.1;
    x2 = -0.1;
    x3 = -0.1;
    
    xesc = zeros(1, N);

    %Start taking time steps
    while any(xesc == 0)
        
        if (n)>=length(nums)
            nums = randn(6, 100000);
            n = 1;
        end

        %Evaluate deterministic part using Heun's scheme.
        kx11 = -(x1^2 - nu1)*(x1 - 1) + betas(1, 2)*(x2-x1) + betas(1, 3)*(x3-x1);
        kx21 = -(x2^2 - nu2)*(x2 - 1) + betas(2, 1)*(x1-x2) + betas(2, 3)*(x3-x2);
        kx31 = -(x3^2 - nu3)*(x3 - 1) + betas(3, 1)*(x1-x3) + betas(3, 2)*(x2-x3);

        %Apply noise
        hkx11 = x1 + h*kx11 + alpha*nums(1, n)*sh;
        hkx21 = x2 + h*kx21 + alpha*nums(2, n)*sh;
        hkx31 = x3 + h*kx31 + alpha*nums(3, n)*sh;

        kx12 = -(hkx11^2 - nu1)*(hkx11 - 1) + betas(1, 2)*(hkx21-hkx11) + betas(1, 3)*(hkx31-hkx11);
        kx22 = -(hkx21^2 - nu2)*(hkx21 - 1) + betas(2, 1)*(hkx11-hkx21) + betas(2, 3)*(hkx31-hkx21);
        kx32 = -(hkx31^2 - nu3)*(hkx31 - 1) + betas(3, 1)*(hkx11-hkx31) + betas(3, 2)*(hkx21-hkx31);

        %Apply heun scheme and noise.
        x1 = x1 + (h/2)*(kx11+kx12) + alpha*nums(4, n)*sh;
        x2 = x2 + (h/2)*(kx21+kx22) + alpha*nums(5, n)*sh;
        x3 = x3 + (h/2)*(kx31+kx32) + alpha*nums(6, n)*sh;

        %Step forward
        tx = tx + h;
        
        if x1 > thresh && xesc(1) == 0
            Exn(1,kn) = tx;
            xesc(1) = 1;
        end
        
        if x2 > thresh && xesc(2) == 0
            Exn(1,kn) = tx;
            xesc(2) = 1;
        end
        
        if x3 > thresh && xesc(3) == 0
            Exn(1,kn) = tx;
            xesc(3) = 1;
        end
        
        if ~any(xesc == 0)
            break;
        end
        
        n = n + 1;

    end
end

%Take average escape time across all realisations for each node
for i = 1:N
    Ex(i) = mean(Exn(:, i));
end

toc