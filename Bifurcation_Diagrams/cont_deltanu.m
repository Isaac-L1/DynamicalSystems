clear
%% Equilibrium Guesses
nu1 = 0.006;
deltanu = 0.001;
nu2 = nu1 + deltanu;
sq1 = sqrt(nu1);
sq2 = sqrt(nu2);

EQ = [-sq1 -sq1 -sq1  sq1  sq1  sq1  1.0  1.0  1.0;
      -sq2  sq2  1.0 -sq2  sq2  1.0 -sq2  sq2  1.0];
  
%% Continuation of equilibrium points for fixed nu

figure(1); clf; hold on;
figure(3); clf; hold on;
for i = 1:9 % Run a continuation for each of the equilibrium points and plot on one figure.
    
    x0 = EQ(:,i);
    pnames = {'nu1', 'deltanu', 'beta'}; %Introduce parameters.
    p0 = [nu1 deltanu 0]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    prob = coco_set(prob, 'cont', 'PtMX', 1000, 'h_max', 0.01, 'h_min', 0.00001);
    ode_args = {@bistuni, x0, pnames, p0};
    cont_args = {1, 'beta', [0 1.35]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation

    figure(1);
    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'beta', 'x', 2); %We continued beta so we need to plot beta against x
    grid on
    figure(3);
    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'beta', 'x'); %We continued beta so we need to plot beta against x
    grid on
end

%% Continuation of the Saddle-Node bifurcation alond beta and nu
% clf(2);
x0 = EQ(:, 2); %Run where we get to the Saddle-Node bifurcation.
pnames = {'nu1', 'deltanu', 'beta'};
p0 = [nu1 deltanu 0];

prob = coco_prob(); %Start a new coco problem
prob = coco_set(prob, 'ode', 'vectorized', false);
%ode_fcns = {@bistuni};
ode_args = {@bistuni, x0, pnames, p0};
cont_args = {1, 'beta', [0 0.35]};

bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation to find the Saddle-Node.

% figure(1); clf; %Show us the specific run that gives us this node.
% thm = struct('special', {{'SN'}});
% coco_plot_bd(thm, 'test2d', 'beta', 'x');
% grid on
    
prob = coco_prob(); %Start a new coco problem
labs = coco_bd_labs(bd1, 'SN'); %Find the saddle node
prob = ode_SN2SN(prob, 'sn1', 'test2d', '', labs(1)); %Continue finding saddle nodes varying nu and beta
bd2 = coco(prob, 'sn', [], 1, {'beta', 'deltanu'}, {[0 1], [-0.01 0.05]});

figure(2); hold on;%Plot the result
thm = struct('special', {{'BP'}});
coco_plot_bd(thm, 'sn', 'beta', 'deltanu');
grid on

%% Continuation of equilibrium points for fixed beta

figure(1); clf; hold on;
for i = 1 % Run a continuation for each of the equilibrium points and plot on one figure.
    
    x0 = EQ(:,i);
    pnames = {'nu1', 'deltanu', 'beta'}; %Introduce parameters.
    p0 = [nu1 deltanu 0.1]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    ode_fcns = {@bistuni};
    ode_args = {ode_fcns{:}, x0, pnames, p0};
    cont_args = {1, 'deltanu', [-0.01 0.05]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation

    figure(1);
    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'x', 'deltanu'); %We continued beta so we need to plot beta against x
    grid on
end

%% Continuation of the Saddle-Node bifurcation alond beta and nu
clf(2);
x0 = EQ(:, 5); %Run where we get to the Saddle-Node bifurcation.
pnames = {'nu1', 'deltanu', 'beta'};
p0 = [nu1 deltanu 0.1];

prob = coco_prob(); %Start a new coco problem
prob = coco_set(prob, 'ode', 'vectorized', false);
ode_fcns = {@bistuni};
ode_args = {ode_fcns{:}, x0, pnames, p0};
cont_args = {1, 'deltanu', [-0.01 0.05]};

bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation to find the Saddle-Node.

% figure(1); clf; %Show us the specific run that gives us this node.
% thm = struct('special', {{'SN'}});
% coco_plot_bd(thm, 'test2d', 'beta', 'x');
% grid on
    
prob = coco_prob(); %Start a new coco problem
labs = coco_bd_labs(bd1, 'SN'); %Find the saddle node
prob = ode_SN2SN(prob, 'sn1', 'test2d', '', labs(2)); %Continue finding saddle nodes varying nu and beta
bd2 = coco(prob, 'sn', [], 1, {'beta', 'deltanu'}, {[0 0.35], [-0.01 0.05]});

figure(2); hold on;%Plot the result
thm = struct('special', {{'BP'}});
coco_plot_bd(thm, 'sn', 'beta', 'deltanu');
grid on



%% Function Declarations

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



