% Critical regions for beta [0 0.006], [0.006 0.066], [0.066 0.204], [0.204
% 0.286], [0.286 0.314], [0.314 Inf)
% Chosen test values beta [0.003 0.03 0.1 0.25 0.3]

beta = 0.3;
alpha = 0.03;
kmax = 2000;
thresh = 0.8;
N = 2;
h = 1e-3;
nus = [0.006; 0.011];
x0 = -sqrt(nus);
Exns = zeros(N, kmax, length(beta));

for b = 1:length(beta)
    betas = [0 beta(b); 0 0];
    Exns(:,:,b) = getEsc(N, nus, betas, h, alpha, x0, kmax, thresh);
end