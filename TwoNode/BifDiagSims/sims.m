% Critical regions for beta [0 0.006], [0.006 0.066], [0.066 0.204], [0.204
% 0.286], [0.286 0.314], [0.314 Inf)
% Chosen test values beta [0.003 0.03 0.1 0.25 0.3]
% dn 0.0001 betas [0.003 0.08 0.15 0.17 0.25 0.4]

beta = linspace(0,0.4,40);
alpha = 0.03;
kmax = 5000;
thresh = 0.8;
N = 2;
h = 1e-3;
nw = 64;
nus = [0.006; 0.0061];
x0 = -sqrt(nus);
Exns1 = zeros(N, kmax, length(beta));
Exns2 = zeros(N, kmax, length(beta));

for b = 1:length(beta)
    betas = [0 beta(b); 0 0];
    p0 = {nus, betas, sum(betas, 2)};
    Exns1(:,:,b) = parGetEsc(@bistable, x0, p0, kmax, thresh, alpha, h, nw);
end

csvwrite('bifdiag1x.csv', reshape(Exns1(1, :,:), [kmax, length(beta)]));
csvwrite('bifdiag1y.csv', reshape(Exns1(2, :,:), [kmax, length(beta)]));

nus = [0.006; 0.011];
x0 = -sqrt(nus);

for b = 1:length(beta)
    betas = [0 beta(b); 0 0];
    p0 = {nus, betas, sum(betas, 2)};
    Exns2(:,:,b) = parGetEsc(@bistable, x0, p0, kmax, thresh, alpha, h, nw);
end

csvwrite('bifdiag2x.csv', reshape(Exns2(1, :,:), [kmax, length(beta)]));
csvwrite('bifdiag2y.csv', reshape(Exns2(2, :,:), [kmax, length(beta)]));

while 1
    pause(1);
end

function f = bistable(x, p)

f = -(x - 1).*(x.^2 - p{1}) + p{2}*x - p{3}.*x;

end