betas = [0, 0, 0; 0, 0, 0; 0, 0, 0];
alpha = 0.03;
nus = [0.01; 0.01; 0.01];
N = 3;
h = 1e-3;
x0 = [-0.1; -0.1; -0.1];
kmax = 2000;
thresh = 0.8;

%Set up escape time matrices and maximum time.
Ex = zeros(N, 1);
Exn = zeros(N, kmax);
B = sum(betas, 2);

%Iterate through number of realisations.
for kn = 1:kmax
    nums = randn(N, 1000000);
    %Set up vectors.
    tx = 0;
    n = 1;
    xesc = zeros(N, 1);
    xn = x0;

    %Start taking time steps
    while any(xesc == 0)
        
        if (2*n > length(nums))
            nums = randn(N, 100000);
            n = 1;
        end
        
        %Find the k1 vector
        k1 = -(xn - 1).*(xn.^2 - nus) + betas*xn - B.*xn;

        %Incorporate noise
        xm = xn + h*k1 + alpha*sqrt(h)*nums(:, 2*n-1);

        %Find the k2 vector
        k2 = -(xm - 1).*(xm.^2 - nus) + betas*xm - B.*xm;

        %Step forward plus noise
        xn = xn + (h/2)*(k1+k2) + alpha*sqrt(h)*nums(:, 2*n);

        %Step forward time and step number
        tx = tx + h;
        n = n + 1;

        ind = (xn > thresh) & (xesc==0);
        xesc(ind) = 1;
        Exn(ind, kn) = tx;


    end
    
end

csvwrite(['sims' simnumber '.csv'], Exn);
