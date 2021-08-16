pl = tiledlayout(2,2);

%% Equilibrium Guesses
nu1 = 0.006;
deltanu = 0.005;
beta = 0.23;
nu2 = nu1 + deltanu;
sq1 = sqrt(nu1);
sq2 = sqrt(nu2);

EQ = [-sq1 -sq1 -sq1  sq1  sq1  sq1  1.0  1.0  1.0;
      -sq2  sq2  1.0 -sq2  sq2  1.0 -sq2  sq2  1.0];
  
%%
nexttile; hold on;
% Continuation of equilibrium points for fixed nu

for i = [1 2 3 7 8 9] % Run a continuation for each of the equilibrium points and plot on one figure.
    
    x0 = EQ(:,i);
    pnames = {'nu1', 'deltanu', 'beta'}; %Introduce parameters.
    p0 = [nu1 deltanu 0]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    prob = coco_set(prob, 'cont', 'PtMX', 1000, 'h_max', 0.01, 'h_min', 0.00001);
    ode_args = {@bistuni, x0, pnames, p0};
    cont_args = {1, 'beta', [0 0.5]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation

    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'beta', 'x'); %We continued beta so we need to plot beta against x
    grid on
end
plot(zeros(1, 13)+beta, (-2:10)./10);

nexttile; hold on;

for i = [1 2 3 7 8 9] % Run a continuation for each of the equilibrium points and plot on one figure.
    
    x0 = EQ(:,i);
    pnames = {'nu1', 'deltanu', 'beta'}; %Introduce parameters.
    p0 = [nu1 deltanu 0]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    prob = coco_set(prob, 'cont', 'PtMX', 1000, 'h_max', 0.01, 'h_min', 0.00001);
    ode_args = {@bistuni, x0, pnames, p0};
    cont_args = {1, 'beta', [0 0.5]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation

    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'beta', 'x', 2); %We continued beta so we need to plot beta against x
    grid on
end
plot(zeros(1, 13)+beta, (-2:10)./10);

%% Continuation of equilibrium points for fixed beta

nexttile; hold on;
for i = [1 2 4 5 8 9] % Run a continuation for each of the equilibrium points and plot on one figure.
    
    x0 = EQ(:,i);
    pnames = {'nu1', 'deltanu', 'beta'}; %Introduce parameters.
    p0 = [nu1 0 beta]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    ode_fcns = {@bistuni};
    ode_args = {ode_fcns{:}, x0, pnames, p0};
    cont_args = {1, 'deltanu', [-0.01 0.02]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation
    
    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'deltanu', 'x'); %We continued beta so we need to plot beta against x
    grid on
end

nexttile; hold on;
for i = [1 2 4 5 8 9] % Run a continuation for each of the equilibrium points and plot on one figure.
    
    x0 = EQ(:,i);
    pnames = {'nu1', 'deltanu', 'beta'}; %Introduce parameters.
    p0 = [nu1 0 beta]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    ode_fcns = {@bistuni};
    ode_args = {ode_fcns{:}, x0, pnames, p0};
    cont_args = {1, 'deltanu', [-0.01 0.02]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation
    
    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'deltanu', 'x', 2); %We continued beta so we need to plot beta against x
    grid on
end




function f = bistbi(x, p)

nu1 = p(1,:);
nu2 = nu1 + p(2,:);
beta = p(3,:);


x1 = x(1,:);
x2 = x(2,:);

f = [(-(x1 - 1.0).*(x1.^2 - nu1) + beta.*(x2-x1));
    (-(x2 - 1.0).*(x2.^2 - nu2) + beta.*(x1-x2))];

end

function f = bistuni(x, p)

nu1 = p(1,:);
nu2 = nu1 + p(2,:);
beta = p(3,:);


x1 = x(1,:);
x2 = x(2,:);

f = [(-(x1 - 1.0).*(x1.^2 - nu1) + beta.*(x2-x1));
    (-(x2 - 1.0).*(x2.^2 - nu2))];

end