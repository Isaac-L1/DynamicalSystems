function roots = cocoEqs(func, par, range, val, tests, pars, ivals)
% COCOEQS - A function to use the coco bifurcation diagram to find roots in
% a given range, given you provide initial guesses sufficient to reach all
% brances. Parameters are:
% - func - function to be zeroed
% - par - parameter to be tested, must be a member of pars
% - range - range to test the parameter in [a b]
% - val - value of par that we want to find equilibria for
% - tests - initial values [w, x; y, z]
% - pars - cell array of parameters {'a', 'b'}
% - ivals - initial values of pars [a0 b0]


roots = [];
n = 1;

for i = 1:length(tests) % Run a continuation for each of the equilibrium points and plot on one figure.

    x0 = tests(:,i);

    prob = coco_prob(); %Open a coco problem
    prob = coco_set(prob, 'ode', 'vectorized', false);
    prob = coco_set(prob, 'cont', 'PtMX', 1000, 'h_max', 0.01, 'h_min', 0.00001);
    prob = ode_isol2ep(prob, '', func, x0, pars, ivals);
    prob = coco_add_event(prob, 'EQ', par, val);

    bd1 = coco(prob, 'test2d', '', 1, par, range); %Run the continuation
    
    labs = coco_bd_labs(bd1, 'EQ');
    
    for l = labs
        roots(1:2, n) = coco_bd_val(bd1, l, 'x');
        eigs = coco_bd_val(bd1, l, 'eigs');
        if eigs(1) < 0 && eigs(2) < 0
            roots(3, n) = 0;
        elseif eigs(1) > 0 && eigs(2) > 0
            roots(3, n) = 2;
        else
            roots(3, n) = 1;
        end
        n = n + 1;
    end

end

end