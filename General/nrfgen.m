function [roots] = nrfgen(N, beta, nu, trials)
%nrfgen - General function to use the higher dimension Newton-Raphson
%method to find the roots of a dynamical system using our model.
% Parameters are:
% - N - number of nodes
% - beta - beta value to test
% - nu - nu value to test
% - trials - initial starting points

f = cell(1, N);
fxn = cell(N);

for i = 1:N
    f{i} = @(xi, x, beta, nu) -xi^3 + xi^2 + xi*nu - nu + beta*sum(x-xi);
end

for i = 1:N
    for j = 1:N
        if i == j
            fxn{i, j} = @(xi, nu, beta) -3.*xi.^2 + 2.*xi + nu - (N-1)*beta;
        else
            fxn{i, j} = @(xi, nu, beta) beta;
        end
    end
end

roots = zeros(N, length(trials));

for i = 1:length(trials)
    x0 = trials(:, i);
    err = 1;
    n = 0;
    x1 = zeros(N, 1);
    while err > 0.0001 && n < 100000
        J = zeros(N);
        fx = zeros(N, 1);
        for j = 1:N
            for k = 1:N
                fxx = fxn{j, k};
                J(j, k) = fxx(x0(j), beta, nu);
            end
            fj = f{j};
            fx(j) = fj(x0(j), x0, beta, nu);
        end
        x1 = x0 - J\fx;
        err = norm(x1-x0);
        x0 = x1;
        n = n + 1;
    end
    roots(:, i) = x1;
end

end