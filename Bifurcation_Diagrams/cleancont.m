clear
%% Equilibrium Guesses
EQ = [-0.1 -0.1 -0.1  0.1  0.1  0.1  1.0  1.0  1.0;
      -0.1  0.1  1.0 -0.1  0.1  1.0 -0.1  0.1  1.0];
  
%% Continuation of equilibrium points for fixed nu

figure(1); clf; hold on;
for i = 1:9 % Run a continuation for each of the equilibrium points and plot on one figure.
    
    x0 = EQ(:, i);
    pnames = {'nu', 'beta'}; %Introduce parameters.
    p0 = [0.01 0]; %Initial conditions on the parameters.

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    ode_fcns = {@bist};
    ode_args = {ode_fcns{:}, x0, pnames, p0};
    cont_args = {1, 'beta', [0 0.25]}; %We want to continue beta and leave nu inactive

    bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation

    figure(1);
    thm = struct('special', {{'SN'}});
    coco_plot_bd(thm, 'test2d', 'beta', 'x'); %We continued beta so we need to plot beta against x
    grid on
end

%% Continuation of the Saddle-Node bifurcation alond beta and nu

x0 = EQ(:, 3); %Run where we get to the Saddle-Node bifurcation.
pnames = {'nu', 'beta'};
p0 = [0.01 0];

prob = coco_prob(); %Start a new coco problem
prob = coco_set(prob, 'ode', 'vectorized', false);
ode_fcns = {@bist};
ode_args = {ode_fcns{:}, x0, pnames, p0};
cont_args = {1, 'beta', [0 0.5]};

bd1 = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:}); %Run the continuation to find the Saddle-Node.

figure(3); clf; %Show us the specific run that gives us this node.
thm = struct('special', {{'SN'}});
coco_plot_bd(thm, 'test2d', 'beta', 'x');
grid on
    
prob = coco_prob(); %Start a new coco problem
labs = coco_bd_labs(bd1, 'SN'); %Find the saddle node
prob = ode_SN2SN(prob, 'sn1', 'test2d', '', labs(1)); %Continue finding saddle nodes varying nu and beta
bd2 = coco(prob, 'sn', [], 1, {'beta', 'nu'}, {[0 1], [0 1]});

figure(2); hold on; %Plot the result
thm = struct('special', {{'BP'}});
coco_plot_bd(thm, 'sn', 'beta', 'nu');
plot(linspace(0,0.25,100), zeros(1,100)+0.01)
grid on

%% Function Declarations

function f = bist(x, p)

nu1 = p(1,:);
nu2 = nu1 + 1e-7;
beta = p(2,:);

x1 = x(1,:);
x2 = x(2,:);

f = [(-(x1 - 1.0).*(x1.^2 - nu1) + beta.*(x2-x1));
    (-(x2 - 1.0).*(x2.^2 - nu2) + beta.*(x1-x2))];

end




