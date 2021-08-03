function [Ex, xr] = getEsc(N, nus, betas, h, alpha, x0, kmax, thresh)
%getEsc - A function which returns an escape time for the parameters
%passed in. 
% Parameters are: 
% - N - number of nodes
% - nus - nu values for each node
% - betas - topology of network with corresponding coupling strength
% - h - step size
% - alpha - amplitude of white noise
% - x0 - vector of initial conditions
% - kmax - number of realisations to use
% - thresh - threshold to use for escape time.
    
    %Check that inputs are in the right form. If not correct.
    if isrow(x0)
        x0 = transpose(x0);
    end
    
    if isrow(nus)
        nus = transpose(nus);
    end
    
    if all(size(nus) ~= [N, 1]) || all(size(betas) ~= [N, N]) || all(size(x0) ~= [N, 1])
        fprintf('Error occured. Check all matrices are of the correct dimensions.\n');
        Ex = zeros(1, N);
        return;
    end
    
    %Set up escape time matrices and maximum time.
    Ex = zeros(N, 1);
    Exn = zeros(N, kmax);
    B = sum(betas, 2);
    
    rng(1)
    %Iterate through number of realisations.
    for kn = 1:kmax
        tic
        %Set up vectors.
        T = 100;
        x = zeros(N, T/h);
        x(:, 1) = x0;
        tx = 0;
        n = 1;
        xesc = zeros(1, N);
        
        %Start taking time steps
        while 1
            
            xn = x(:, n);
        
            %Find the k1 vector
            k1 = -(xn - 1).*(xn.^2 - nus) + betas*xn - B.*xn;

            %Incorporate noise
            xm = xn + h*k1 + alpha*sqrt(h)*randn(N, 1);

            %Find the k2 vector
            k2 = -(xm - 1).*(xm.^2 - nus) + betas*xm - B.*xm;

            %Step forward plus noise
            x(:, n+1) = xn + (h/2)*(k1+k2) + alpha*sqrt(h)*randn(N, 1);

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
        toc
    end
    
    %Take average escape time across all realisations for each node
    for i = 1:N
        Ex(i) = mean(Exn(:, i));
    end

end