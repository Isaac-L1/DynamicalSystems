beta = 0;
nu1 = 0.01;
nu2 = 0.01;

f = @(x, y, nu, beta) -x^3 + x^2 + x*nu - nu + beta*y - beta*x;
g = @(x, y, nu, beta) -y^3 + y^2 + y*nu - nu + beta*x - beta*y;
fx = @(x, nu, beta) -3*x^2 + 2*x + nu - beta;
fy = beta;
gx = beta;
gy = @(y, nu, beta) -3*y^2 + 2*y + nu - beta;

trials = [-1, -1, -1, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5; -1, 0.5, 1.5, -1, 0.5, 1.5, -1, 0.5, 1.5];
%trials = [0; 0];
roots = zeros(3, 9);
eigs = zeros(3, 2, 9);

for i = 1:length(trials)
    
    x0 = trials(:, i);
    err = 1;
    n = 0;
    x1 = [0; 0];
    
    while err > 0.0001 && n < 100000
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