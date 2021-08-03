betas = [0, 0.1; 0, 0];
alpha = 0.03;
meanNu = 0.0001;
dN1 = 0.007;
dN2 = 0.006;
err3 = 1;
firstRun = true;
secondRun = true;
deltaNu = dN1;
err1 = 1;
err2 = 1;

tic

while abs(err3) > 0.001
    
    nus = [meanNu; meanNu + deltaNu];
    N = 2;
    h = 1e-3;
    x0 = [-sqrt(nus(1)); -sqrt(nus(2))];
    kmax = 10000;
    thresh = 0.8;

    %Set up escape time matrices and maximum time.
    Exn = zeros(N, kmax);
    B = sum(betas, 2);

    %Iterate through number of realisations.
    parfor kn = 1:kmax
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

            for i = 1:N
                if xn(i) > thresh && xesc(i) == 0
                    xesc(i) = 1;
                    Exn(i, kn) = tx;
                end
            end
        end
    end
    
    Emin = min(Exn);
    P = sum(Exn(2, :) == Emin)/kmax;
    err3 = P-0.5;
    
    if firstRun
        deltaNu = dN2;
        firstRun = false;
        err1 = err3;
        continue;
    end
    
    if secondRun
        secondRun = false;
        if err1*err3>0
            error("Suggested means do not bracket P=0.5.");
        end
        deltaNu = (dN1+dN2)/2;
        err2 = err3;
        continue;
    end
    
    if err1*err3 < 0
        dN2 = deltaNu;
        deltaNu = (dN1+deltaNu)/2;
        err2 = err3;
    elseif err2*err3 < 0
        dN1 = deltaNu;
        deltaNu = (dN2+deltaNu)/2;
        err1 = err3;
    else
        error("An error occured while evaluating which side of the bisection to follow.");
    end
    
    fprintf(['\n Now testing a deltaNu = ' num2str(deltaNu) '. \n']);
    
end

finaldN = deltaNu;

fprintf(['\n Root found at deltaNu = ' num2str(deltaNu) '. \n']);

toc





