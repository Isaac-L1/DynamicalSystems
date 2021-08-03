function [roots, eigs] = nrfuncuni(beta1, beta2, nu1, nu2)
%nrfunc - Function to use the 2D Newton-Raphson method to find the roots of
%the model for the two node bidirectionally coupled case.
% The parameters are:
% - beta - beta value to test
% - nu - nu value to test

%Problem specific functions
f = @(x, y, nu, beta) -x^3 + x^2 + x*nu - nu + beta*y - beta*x;
g = @(x, y, nu, beta) -y^3 + y^2 + y*nu - nu + beta*x - beta*y;
fx = @(x, nu, beta) -3*x^2 + 2*x + nu - beta;
fy = beta1;
gx = beta2;
gy = @(y, nu, beta) -3*y^2 + 2*y + nu - beta;

sq1 = sqrt(nu1);
sq2 = sqrt(nu2);

trials = linspace(-0.5, 1.5, 20);
%trials = [0; 0];
roots = zeros(3, length(trials), length(trials));
eigs = zeros(3, 2, length(trials), length(trials));

%One run for each trial
for i = 1:length(trials)
    for j = 1:length(trials)
    
        x0 = [trials(i); trials(j)];
        err = 1;
        n = 0;
        x1 = [0; 0];

        %Iterate through the newton raphson method
        while err > 1e-6
            J = [fx(x0(1), nu1, beta1), fy; gx, gy(x0(2), nu2, beta2)];
            x1 = x0 - J\[f(x0(1), x0(2), nu1, beta1); g(x0(1), x0(2), nu2, beta2)];
            err = norm(x1-x0);
            x0 = x1;
            n = n + 1;
        end

        %Characterise the equilibrium
        roots(1:2, i, j) = x1;
        [vects, vals] = eig(J);
        if vals(1, 1)<0 && vals(2, 2)<0
            roots(3, i, j) = 1;
            eigs(3, :, i, j) = 1;
        elseif vals(1, 1)>0 && vals(2, 2)>0
            roots(3, i, j) = 2;
            eigs(3, :, i, j) = 2;
        elseif vals(1, 1)>0
            roots(3, i, j) = 3;
            eigs(3, 1, i, j) = 2;
            eigs(3, 2, i, j) = 1;
        else
            roots(3, i, j) = 3;
            eigs(3, 1, i, j) = 1;
            eigs(3, 2, i, j) = 2;
        end

        eigs(1:2, 1, i, j) = vects(:, 1);
        eigs(1:2, 2, i, j) = vects(:, 2);
    
    end
    
end

end