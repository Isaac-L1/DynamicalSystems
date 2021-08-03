f = @(x, y, nu, beta) -x^3 + x^2 + x*nu - nu + beta*y - beta*x;
g = @(y, nu) -y^3 + y^2 + y*nu - nu;
fx = @(x, nu, beta) -3*x^2 + 2*x + nu - beta;
fy = 0.1;
gx = 0;
gy = @(y, nu, beta) -3*y^2 + 2*y + nu - beta;

% Find y
nu1 = 0.005;
nu2 = 0.007;

y1 = 0.08;
y2 = 0.09;
r1 = g(y1,nu2);
r2 = g(y2,nu2);
err = 1;

while abs(err) > 1e-8
    yt = (y1+y2)/2;
    rt = g(yt,nu2);
    err = rt;
    if r1 * rt < 0
        y2 = yt;
        r2 = rt;
    elseif r2 * rt < 0
        y1 = yt;
        r1 = rt;
    else
        error('You plonker')
    end
end

x1 = -0.05;
x2 = 0.2;
beta = 0.01;
r1 = f(x1, yt, nu1, beta);
r2 = f(x2, yt, nu1, beta);
err = 1;

while abs(err) > 1e-8
    xt = (x1+x2)/2;
    rt = f(xt, yt, nu1, beta);
    err = rt;
    if r1 * rt < 0
        x2 = xt;
        r2 = rt;
    elseif r2 * rt < 0
        x1 = xt;
        r1 = rt;
    else
        error('You plonker')
    end
end


