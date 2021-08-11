% dn 0.005 betas [0.003 0.03 0.1 0.25 0.3 0.5]
% dn 0.0001 betas [0.003 0.08 0.15 0.17 0.25 0.4]

pl = tiledlayout(2,2);
pln = 1;

nu1 = 0.006;
deltanu = 0.005;
nu2 = nu1 + deltanu;
sq1 = sqrt(nu1);
sq2 = sqrt(nu2);

EQ = [-sq1 -sq1 -sq1  sq1  sq1  sq1  1.0  1.0  1.0;
      -sq2  sq2  1.0 -sq2  sq2  1.0 -sq2  sq2  1.0];

nexttile(pln); hold on;
for i = 1:9 % Run a continuation for each of the equilibrium points and plot on one figure.

    x0 = EQ(:,i);
    pnames = {'nu1', 'deltanu', 'beta'}; %Introduce parameters.
    p0 = [nu1 deltanu 0]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    prob = coco_set(prob, 'cont', 'PtMX', 1000, 'h_max', 0.01, 'h_min', 0.00001);
    ode_args = {@bistuni, x0, pnames, p0};
    cont_args = {1, 'beta', [0 0.4]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation

    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'beta', 'x'); %We continued beta so we need to plot beta against x
    grid on
end
plot(zeros(1,12)+0.003, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.03, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.1, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.25, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.3, -0.1:0.1:1, 'color', [0 0 0 0.5]);
title(['Bifurcation Diagram for \delta \nu = ' num2str(deltanu)]);
xlabel('\beta');
pln = pln + 1;

%%
nu1 = 0.006;
deltanu = 0.0001;
nu2 = nu1 + deltanu;
sq1 = sqrt(nu1);
sq2 = sqrt(nu2);

EQ = [-sq1 -sq1 -sq1  sq1  sq1  sq1  1.0  1.0  1.0;
      -sq2  sq2  1.0 -sq2  sq2  1.0 -sq2  sq2  1.0];

nexttile(pln); hold on;
for i = 1:9 % Run a continuation for each of the equilibrium points and plot on one figure.

    x0 = EQ(:,i);
    pnames = {'nu1', 'deltanu', 'beta'}; %Introduce parameters.
    p0 = [nu1 deltanu 0]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    prob = coco_set(prob, 'cont', 'PtMX', 1000, 'h_max', 0.01, 'h_min', 0.00001);
    ode_args = {@bistuni, x0, pnames, p0};
    cont_args = {1, 'beta', [0 0.4]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation

    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'beta', 'x'); %We continued beta so we need to plot beta against x
    grid on
end
plot(zeros(1,12)+0.003, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.08, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.15, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.17, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.25, -0.1:0.1:1, 'color', [0 0 0 0.5]);
plot(zeros(1,12)+0.4, -0.1:0.1:1, 'color', [0 0 0 0.5]);
title(['Bifurcation Diagram for \delta \nu = ' num2str(deltanu)]);
xlabel('\beta');
pln = pln + 1;

%%
load('simdata1.mat');
EDS = zeros(2, 5);
betas = [0.003 0.03 0.1 0.25 0.3];

for j = 1:5
    Exn = Exns(:,:,j);

    X1 = Exn(1,:);
    X2 = Exn(2,:);
    Emin = min(Exn);
    Emax = max(Exn);

    Ex1 = mean(Exn(1, :));
    Ex2 = mean(Exn(2, :));

    Sx1 = std(Exn(1, :));
    Sx2 = std(Exn(2, :));

    pos = zeros(2, length(Exn));
    type = zeros(1, length(Exn));
    ED = zeros(2, length(Exn));

    for i = 1:length(Exn)
        mn = find(Exn(:, i) == Emin(i));
        pos(mn, i) = 1;
        mx = find(Exn(:, i) == Emax(i));
        pos(mx, i) = 2;

    end

    EO = sort(Exn);

    ED(1, :) = EO(1, :);
    ED(2, :) = EO(2, :) - EO(1, :);

    indX = (pos(1, :) == 1);
    PX = sum(indX(:)==1)/10000;
    EDX = ED(:, indX);
    indY = (pos(2, :) == 1);
    PY = sum(indY(:)==1)/10000;
    EDY = ED(:, indY);
    EED = mean(ED, 2);
    SED = std(ED, 0, 2);
    EEDX = mean(EDX, 2);
    SEDX = std(EDX, 0, 2);
    EEDY = mean(EDY, 2);
    SEDY = std(EDY, 0, 2);

    CED = SED./EED;
    CEDX = SEDX./EEDX;
    CEDY = SEDY./EEDY;
    
    EDS(:,j) = EED;
end

nexttile(pln); hold on;
plot(betas, EDS(1,:));
plot(betas, EDS(2,:));
plot(0:0.1:0.4, zeros(1,5), 'k');
title('Mean first and domino escape times for tested critical regions');
xlabel('\beta');
ylabel('E(\tau)');
pln = pln + 1;

%%
load('simdata12.mat');
EDS = zeros(2, 20);
betas = linspace(0,0.4,20);

for j = 1:20
    Exn = Exns(:,:,j);

    X1 = Exn(1,:);
    X2 = Exn(2,:);
    Emin = min(Exn);
    Emax = max(Exn);

    Ex1 = mean(Exn(1, :));
    Ex2 = mean(Exn(2, :));

    Sx1 = std(Exn(1, :));
    Sx2 = std(Exn(2, :));

    pos = zeros(2, length(Exn));
    type = zeros(1, length(Exn));
    ED = zeros(2, length(Exn));

    for i = 1:length(Exn)
        mn = find(Exn(:, i) == Emin(i));
        pos(mn, i) = 1;
        mx = find(Exn(:, i) == Emax(i));
        pos(mx, i) = 2;

    end

    EO = sort(Exn);

    ED(1, :) = EO(1, :);
    ED(2, :) = EO(2, :) - EO(1, :);

    indX = (pos(1, :) == 1);
    PX = sum(indX(:)==1)/10000;
    EDX = ED(:, indX);
    indY = (pos(2, :) == 1);
    PY = sum(indY(:)==1)/10000;
    EDY = ED(:, indY);
    EED = mean(ED, 2);
    SED = std(ED, 0, 2);
    EEDX = mean(EDX, 2);
    SEDX = std(EDX, 0, 2);
    EEDY = mean(EDY, 2);
    SEDY = std(EDY, 0, 2);

    CED = SED./EED;
    CEDX = SEDX./EEDX;
    CEDY = SEDY./EEDY;
    
    EDS(:,j) = EED;
end

nexttile(pln); hold on;
plot(betas, EDS(1,:));
plot(betas, EDS(2,:));
plot(0:0.1:0.4, zeros(1,5), 'k');
title('Mean first and domino escape times for tested critical regions');
xlabel('\beta');
ylabel('E(\tau)');
%%

function f = bistuni(x, p)

nu1 = p(1,:);
nu2 = nu1 + p(2,:);
beta = p(3,:);


x1 = x(1,:);
x2 = x(2,:);

f = [(-(x1 - 1.0).*(x1.^2 - nu1) + beta.*(x2-x1));
    (-(x2 - 1.0).*(x2.^2 - nu2))];

end