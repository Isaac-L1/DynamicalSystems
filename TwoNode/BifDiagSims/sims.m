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
nus = [0.006; 0.0061];
x0 = -sqrt(nus);
Exns1 = zeros(N, kmax, length(beta));

parfor b = 1:length(beta)
    betas = [0 beta(b); 0 0];
    Exns1(:,:,b) = getEsc(N, nus, betas, h, alpha, x0, kmax, thresh);
end

beta = linspace(0,0.4,40);
alpha = 0.03;
kmax = 10000;
thresh = 0.8;
N = 2;
h = 1e-3;
nus = [0.006; 0.011];
x0 = -sqrt(nus);
Exns2 = zeros(N, kmax, length(beta));

parfor b = 1:length(beta)
    betas = [0 beta(b); 0 0];
    Exns2(:,:,b) = getEsc(N, nus, betas, h, alpha, x0, kmax, thresh);
end

while 1
    pause(1);
end