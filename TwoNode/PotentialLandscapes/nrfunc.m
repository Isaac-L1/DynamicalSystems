function [roots, eigs] = nrfunc(beta, nu1, nu2)
%nrfunc - Function to use the 2D Newton-Raphson method to find the roots of
%the model for the two node bidirectionally coupled case.
% The parameters are:
% - beta - beta value to test
% - nu - nu value to test

f = @(x, y, nu, beta) -x^3 + x^2 + x*nu - nu + beta*y - beta*x;
g = @(x, y, nu, beta) -y^3 + y^2 + y*nu - nu + beta*x - beta*y;
fx = @(x, nu, beta) -3*x^2 + 2*x + nu - beta;
fy = beta;
gx = beta;
gy = @(y, nu, beta) -3*y^2 + 2*y + nu - beta;

sq1 = sqrt(nu1);
sq2 = sqrt(nu2);

trials = [-sq1 -sq1 -sq1  sq1  sq1  sq1  1.0  1.0  1.0;
          -sq2  sq2  1.0 -sq2  sq2  1.0 -sq2  sq2  1.0];
%trials = [0; 0];
roots = zeros(3, length(trials));
eigs = zeros(3, 2, length(trials));

for i = 1:length(trials)
    
    x0 = trials(:, i);
    err = 1;
    n = 0;
    x1 = [0; 0];
    
    while n < 100
        J = [fx(x0(1), nu1, beta), fy; gx, gy(x0(2), nu2, beta)];
        x1 = x0 - J\[f(x0(1), x0(2), nu1, beta); g(x0(1), x0(2), nu2, beta)];
        err = norm(x1-x0);
        x0 = x1;
        n = n + 1;
    end
    
    roots(1:2, i) = x1;
    [vects, vals] = eig(J);
    if vals(1, 1)<0 && vals(2, 2)<0
        roots(3, i) = 1;
        eigs(3, :, i) = 1;
    elseif vals(1, 1)>0 && vals(2, 2)>0
        roots(3, i) = 2;
        eigs(3, :, i) = 2;
    elseif vals(1, 1)>0
        roots(3, i) = 3;
        eigs(3, 1, i) = 2;
        eigs(3, 2, i) = 1;
    else
        roots(3, i) = 3;
        eigs(3, 1, i) = 1;
        eigs(3, 2, i) = 2;
    end
    
    eigs(1:2, 1, i) = vects(:, 1);
    eigs(1:2, 2, i) = vects(:, 2);
    
end

end