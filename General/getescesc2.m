betasm = [0, 0; 0, 0];
alpha = 0.1;
nus = [0.01; 0.01];
N = 2;
h = 1e-3;
T = 150;
x0 = zeros(2*sims, 1);
kmax = 10;
thresh = 0.8;
sims = 10;
betas = zeros(2*sims);
for i = 1:sims
    betas(2*i-1:2*i,2*i-1:2*i) = betasm;
end

%Set up escape time matrices and maximum time.
Ex = zeros(N, 1);
Exn = zeros(N, kmax);
B = sum(betas, 2);
tic

%Iterate through number of realisations.
for kn = 1:kmax
    nums = randn(N, 2*T/h);
    %Set up vectors.
    x = zeros(N, T/h);
    x(:, 1) = x0;
    tx = 0;
    n = 1;
    xesc = zeros(1, N);
    xn = x0;

    %Start taking time steps
    while 1
        
        %Find the k1 vector
        k1 = -(xn - 1).*(xn.^2 - nus) + betas*xn - B.*xn;

        %Incorporate noise
        xm = xn + h*k1 + alpha*sqrt(h)*nums(:, 2*n-1);

        %Find the k2 vector
        k2 = -(xm - 1).*(xm.^2 - nus) + betas*xm - B.*xm;

        %Step forward plus noise
        xn1 = xn + (h/2)*(k1+k2) + alpha*sqrt(h)*nums(:, 2*n);
        
        x(:, n+1) = xn1;
        xn = xn1;

        %Step forward time and step number
        tx = tx + h;
        n = n + 1;

        %Check which nodes have escaped
%             xesc = floor(x(:, n)./thresh);
%             Exn(:, kn) = tx.*xesc;
        for i = 1:N
            if (x(i, n) > thresh) && xesc(i) == 0
                xesc(i) = 1;
                Exn(kn, i) = tx;
            end
        end

        %If all nodes have escaped then exit for loop
        if all(xesc == 1)
            break
        end

    end
    xr = x;
    
end
toc
%Take average escape time across all realisations for each node
for i = 1:N
    Ex(i) = mean(Exn(:, i));
end