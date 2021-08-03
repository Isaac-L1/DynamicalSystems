clear
%% Equilibrium Guesses

pl = tiledlayout(2, 4);
pln = 1;
dns = [0 1e-4 1e-3 1e-2];
  
%% Continuation of equilibrium points for fixed nu

for j = dns
    nu1 = 0.006;
    deltanu = j;
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
        ode_args = {@bistbi, x0, pnames, p0};
        cont_args = {1, 'beta', [0 0.3]}; %We want to continue beta and leave nu inactive

        bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation
        
        figure(1);
        thm = struct('special', {{'SN'}});
        coco_plot_bd(thm, 'test2d', 'beta', 'x'); %We continued beta so we need to plot beta against x
        grid on
    end
    title(['\delta \nu = ' num2str(deltanu)]);
    pln = pln + 1;
end

for j = dns
    nu1 = 0.006;
    deltanu = j;
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
        ode_args = {@bistbi, x0, pnames, p0};
        cont_args = {1, 'beta', [0 0.3]}; %We want to continue beta and leave nu inactive

        bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation
        
        figure(1);
        thm = struct('special', {{'SN'}});
        coco_plot_bd(thm, 'test2d', 'beta', 'x', 2); %We continued beta so we need to plot beta against x
        grid on
    end
    title(['\delta \nu = ' num2str(deltanu)]);
    pln = pln + 1;
end
title(pl, 'Bifurcation Diagrams for \nu 1 = 0.006');

%%
function f = bistbi(x, p)

nu1 = p(1,:);
nu2 = nu1 + p(2,:);
beta = p(3,:);


x1 = x(1,:);
x2 = x(2,:);

f = [(-(x1 - 1.0).*(x1.^2 - nu1) + beta.*(x2-x1));
    (-(x2 - 1.0).*(x2.^2 - nu2) + beta.*(x1-x2))];

end

