function [Exn] = getEsc2(N, Nus, Betas, h, alpha, x0, kmax, thresh)
%getEsc - A function which returns an escape time for the parameters
%passed in. 
% Parameters are: 
% - N - number of nodes
% - Nus - nu values for each node
% - Betas - topology of network with corresponding coupling strength
% - h - step size
% - alpha - amplitude of white noise
% - X0 - vector of initial conditions
% - kmax - number of realisations to use
% - thresh - threshold to use for escape time.
    
    %Check that inputs are in the right form. If not correct.
    if isrow(x0)
        x0 = transpose(x0);
    end
    
    if isrow(Nus)
        Nus = transpose(Nus);
    end
    
    if all(size(Nus) ~= [N, 1]) || all(size(Betas) ~= [N, N]) || all(size(x0) ~= [N, 1])
        fprintf('Error occured. Check all matrices are of the correct dimensions.\n');
        Exn = zeros(1, N);
        return;
    end
    
    %Set up escape time matrices and maximum time.
    Exn = zeros(N, kmax);
    B = sum(Betas, 2);
    nw = 8;
    dt = h;
    sh = sqrt(h);
    bs = kron(eye(nw), Betas);
    betas = sparse(bs);
    B = repmat(B, [nw,1]);
    nus = repmat(Nus, [nw,1]);
    w = zeros(1, nw);
    nums = randn(N*nw, 100000);
    completed = zeros(1, kmax);
    wN = 1:N*nw;
    xn = zeros(N*nw, 1);
    
    tic
    %Set up vectors.
    tx = 0;
    xesc = zeros(N*nw, 1);
    n = 1;
    knext = 1;

    %Start taking time steps
    while any(completed == 0)
         
        if any(w == 0)
            empty = (w==0);
            ei = find(empty);
            if sum(empty) == 1
                if knext <= kmax
                    w(ei) = knext;
                    xn(N*(ei-1)+1:N*ei) = x0;
                    knext = knext + 1;
                    xesc(N*(ei-1)+1:N*ei) = 0;
                else
                    nw = nw - 1;
                    bs(N*(ei-1)+1:N*ei, :) = [];
                    bs(:, (N*(ei-1)+1:N*ei)) = [];
                    betas = sparse(bs);
                    B(N*(ei-1)+1:N*ei) = [];
                    nus(N*(ei-1)+1:N*ei) = [];
                    w(ei) = [];
                    xesc(N*(ei-1)+1:N*ei) = [];
                    xn(N*(ei-1)+1:N*ei) = [];
                    wN = 1:N*nw;
                    nums = randn(N*nw, 100000);
                    n = 1;
                end
            else
                for e = 1:sum(empty)
                    if knext <= kmax
                        w(ei(e)) = knext;
                        xn(N*(ei(e)-1)+1:N*ei(e)) = x0;
                        knext = knext + 1;
                        xesc(N*(ei(e)-1)+1:N*ei(e)) = 0;
                    else
                        nw = nw - 1;
                        bs(N*(ei(e)-1)+1:N*ei(e), :) = [];
                        bs(:, (N*(ei(e)-1)+1:N*ei(e))) = [];
                        betas = sparse(bs);
                        B(N*(ei(e)-1)+1:N*ei(e)) = [];
                        nus(N*(ei(e)-1)+1:N*ei(e)) = [];
                        w(ei(e)) = [];
                        xesc(N*(ei(e)-1)+1:N*ei(e)) = [];
                        xn(N*(ei(e)-1)+1:N*ei(e)) = [];
                        wN = 1:N*nw;
                        nums = randn(N*nw, 100000);
                        n = 1;
                    end
                end
            end
        end
        
        if (2*n)>=length(nums)
            nums = randn(N*nw, 100000);
            n = 1;
        end

        %Find the k1 vector
        k1 = -(xn - 1).*(xn.^2 - nus) + betas*xn - B.*xn;

        %Incorporate noise
        xm = xn + h.*k1 + alpha.*sh.*nums(:, 2*n - 1);

        %Find the k2 vector
        k2 = -(xm - 1).*(xm.^2 - nus) + betas*xm - B.*xm;

        %Step forward plus noise
        xn = xn + (h/2).*(k1+k2) + alpha.*sh.*nums(:, 2*n);

        %Step forward time and step number
        tx = tx + dt;
        n = n + 1;

        ind = (xn > thresh) & (xesc==0);
        xesc(ind) = 1;
        Exn((mod(wN(ind), N) + 1), w(ceil(wN(ind)./N))) = tx;
        
        xescr = reshape(xesc, [N nw]);
        escind = (~any(xescr == 0));
        completed(w(escind)) = 1;
        w(escind) = 0;

    end
    toc

end