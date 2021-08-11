function Exn = parGetEsc(func, x0, p0, kmax, thresh, alpha, h, nw)
% PARGETESC - Function that finds escape times using a parallelised method,
% generalised for any dynamical system but always implements the Heun
% method. Parameters are:
% - func - function handle to be evaluated. Constraints are listed below.
% - x0 - Nx1 vector of initial x values.
% - p0 - cell array of syatem parameters.
% - kmax - number of realisations to compute.
% - thresh - threshold of x for detecting escape.
% - alpha - noise amplitude
% - h - step size
% - nw - number of workers to paralellise to.
% 
% func must be such that it takes two parameters, x (the node positions) and
% p (a cell array containing system parameters)
% 
% for example:
% 
% function x1 = func(x, p)
% 
% x1 = x.^2 + p{1}.*x + p{2}*x;
% 
% end

    N = length(x0); %Number of nodes (implied from inputs)
    Exn = zeros(N, kmax); %Set up return tensor.
    sh = sqrt(h); %Saves repeated floating point operations in iteration
    w = zeros(1, nw); %Vector which stores which realisation each worker is on
    nums = randn(N*nw, 100000); %Pregenerate random numbers for performance
    completed = zeros(1, kmax); %Store which realisations have completed
    wN = 1:N*nw; %Vector used in calucations in line LINEREF
    
    xn = zeros(N*nw, 1); %Prepare vector to store current position of nodes
    
    tx = zeros(N*nw, 1); %Initial time
    xesc = zeros(N*nw, 1); %Vector which detects node escapes
    n = 1; %Initial step number for random number matrix
    knext = 1; %First realisation to compute
    p = cell(size(p0)); %Prepare cell array to store system parameters
    pd = cell(size(p0)); %Prepare cell array to store dense versions of sparse matrices
    
    for i = 1:length(p0) %For each system parameter, check whether it is a vector or a matrix
        s = size(p0{i});
        if s(2) == 1
            p{i} = repmat(p0{i}, [nw,1]); %If vector simply repeat it into required format
        else
            pd{i} = kron(eye(nw), p0{i}); %If matrix, find kronecker product.
            p{i} = sparse(pd{i}); %Store sparsely
        end
    end

    %Start taking time steps
    while any(completed == 0) %Once all realisations are computed, exit
         
        if any(w == 0) %If any workers are idle, assign realisations
            empty = (w==0);
            ei = find(empty);
            for e = 1:sum(empty)
                if knext <= kmax %If there are realisations left assign.
                    w(ei(e)) = knext;
                    xn(N*(ei(e)-1)+1:N*ei(e)) = x0;
                    tx(N*(ei(e)-1)+1:N*ei(e)) = 0;
                    knext = knext + 1;
                    xesc(N*(ei(e)-1)+1:N*ei(e)) = 0;
                else %Otherwise remove worker from the stack.
                    nw = nw - 1;
                    for i = 1:length(p)
                        s = size(p{i});
                        if s(2) == 1
                            p{i}(N*(ei-1)+1:N*ei) = [];
                        else
                            pd{i}(N*(ei-1)+1:N*ei, :) = [];
                            pd{i}(:, (N*(ei-1)+1:N*ei)) = [];
                            p{i} = sparse(pd{i});
                        end
                    end
                    w(ei(e)) = [];
                    xesc(N*(ei(e)-1)+1:N*ei(e)) = [];
                    xn(N*(ei(e)-1)+1:N*ei(e)) = [];
                    tx(N*(ei(e)-1)+1:N*ei(e)) = [];
                    wN = 1:N*nw;
                    nums = randn(N*nw, 100000);
                    n = 1;
                end
            end
        end
        
        %If we have ran out of random numbers, regenerate.
        if (2*n)>=length(nums)
            nums = randn(N*nw, 100000);
            n = 1;
        end

        %Find the k1 vector
        k1 = func(xn, p);

        %Incorporate noise
        xm = xn + h.*k1 + alpha.*sh.*nums(:, 2*n - 1);

        %Find the k2 vector
        k2 = func(xm, p);

        %Step forward plus noise
        xn = xn + (h/2).*(k1+k2) + alpha.*sh.*nums(:, 2*n);

        %Step forward time and step number
        tx = tx + h;
        n = n + 1;

        %Find nodes that have escaped
        ind = (xn > thresh) & (xesc==0);
        xesc(ind) = 1;
        Exn((mod(wN(ind), N) + 1), w(ceil(wN(ind)./N))) = tx(ind);
        
        %If any whole systems have escaped then make that worker idle.
        xescr = reshape(xesc, [N nw]);
        escind = (~any(xescr == 0));
        completed(w(escind)) = 1;
        w(escind) = 0;

    end

end