betas = [0, 0.1; 0, 0];
alpha = 0.03;
N = 2;
h = 1e-3;
kmax = 10000;
thresh = 0.8;
meanNu = 0.0001;
dN1 = 0.01;
dN2 = 0.004;
err3 = 1;
firstRun = true;
secondRun = true;
deltaNu = dN1;
err1 = 1;
err2 = 1;

tic

while abs(err3) > 0.005
    
    fprintf(['\n Now testing a deltaNu = ' num2str(deltaNu) '. \n']);
    
    nus = [meanNu; meanNu + deltaNu];
    x0 = [-sqrt(nus(1)); -sqrt(nus(2))];

    Exn = getEsc(N, nus, betas, h, alpha, x0, kmax, thresh);
    
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
    
    finaldN = deltaNu;
    
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
    
    
    
end

fprintf(['\n Root found at deltaNu = ' num2str(deltaNu) '. \n']);

toc
