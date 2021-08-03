function [x] = genReal(N, nus, betas, h, T, alpha, x0)
%genReal - A function which returns a realisation vector for the parameters
%passed in.
% Parameters are: 
% - N - number of nodes
% - nus - nu values for each node
% - betas - topology of network with corresponding coupling strength
% - h - step size
% - T - total time to step through
% - alpha - amplitude of white noise
% - x0 - vector of initial conditions

    %Check that inputs are in the right form. If not then correct.
    if isrow(x0)
        x0 = transpose(x0);
    end
    
    if isrow(nus)
        nus = transpose(nus);
    end
    
    if any(size(nus) ~= [N, 1]) || any(size(betas) ~= [N, N]) || any(size(x0) ~= [N, 1])
        fprintf('Error occured. Check all matrices are of the correct dimensions.\n');
        x = zeros(N, T/h+1);
        return;
    end
    
    %Set up Vectors
    x = zeros(N, T/h);
    x(:, 1) = x0;
    tx = 0;
    n = 1;
    
    tic
    %Start making time steps
    while n <= (T/h)
        xn = x(:, n);
        
        %Find the k1 vector
        k1 = -(xn - 1).*(xn.^2 - nus);
        for i = 1:N
            k1(i) = k1(i) + betas(i, :)*(xn-xn(i));
        end
        
        %Incorporate noise
        xm = xn + h*k1 + alpha*sqrt(h)*randn(N, 1);
        
        %Find the k2 vector
        k2 = -(xm - 1).*(xm.^2 - nus);
        for i = 1:N
            k2(i) = k2(i) + betas(i, :)*(xm-xm(i));
        end
        
        %Step forward plus noise
        x(:, n+1) = xn + (h/2)*(k1+k2) + alpha*sqrt(h)*randn(N, 1);
        
        %Step forward time and step number
        tx = tx + h;
        n = n + 1;
    end
    toc

end