N = 2;
alpha = 0.1;
Nus = [0.0001;0.0001];
Betas = [0 0.1; 0.1 0];
x0 = [-0.01; -0.01];
h = 1e-3;
kmax = 10;
thresh = 0.8;

esc = getEsc2(N, Nus, Betas, h, alpha, x0, kmax, thresh);