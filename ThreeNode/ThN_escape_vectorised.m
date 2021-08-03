clear

N = 3;
alpha = 0.1;
betas = [0,0,0;0,0,0;0,0,0];
nus = [0.01;0.01;0.01];
x0 = [-0.1;-0.1;-0.1];
h = 1e-3;
sh = sqrt(h);
kmax = 1000;
thresh = 0.6;

%%
tic

%Set up escape time matrices and maximum time.
Ex = zeros(N, 1);
Exn = zeros(N, kmax);
B = sum(betas, 2);

%Iterate through number of realisations.
for kn = 1:kmax
    
    nums = randn(6, 500000);
    %Set up vectors.
    tx = 0;
    n = 1;
    xn = x0;
    
    xesc = zeros(1, N);

    %Start taking time steps
    while any(xesc == 0)
        
        if (n)>=length(nums)
            nums = randn(6, 100000);
            n = 1;
        end

        %Find the k1 vector
        k1 = -(xn - 1).*(xn.^2 - nus) + betas*xn - B.*xn;

        %Incorporate noise
        xm = xn + h.*k1 + alpha*sh*nums(1:3, n);

        %Find the k2 vector
        k2 = -(xm - 1).*(xm.^2 - nus) + betas*xm - B.*xm;

        %Step forward plus noise
        xn = xn + (h/2).*(k1+k2) + alpha*sh*nums(4:6, n);

        %Step forward time and step number
        tx = tx + h;

        ind = (xn > thresh) & (xesc.'==0);
        xesc(ind) = 1;
        Exn(ind, kn) = tx;
        
        n = n + 1;

    end
end

%Take average escape time across all realisations for each node
for i = 1:N
    Ex(i) = mean(Exn(:, i));
end

toc