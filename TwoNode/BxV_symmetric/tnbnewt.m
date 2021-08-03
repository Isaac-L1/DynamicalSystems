clear

%%
beta = 0;
%nus = logspace(-4, -1.3, 200);
nus1 = 0.001;
nus2 = 0.01;
%noroots = zeros(length(beta), length(nus));
%allroots = zeros(2, 9, length(betas));

tic
for i = 1:length(nus1)
    %fprintf(['\n Now calculating roots for beta = ' num2str(betas(i)) '. \n'])
    for j = 1:length(nus2)
        roots = nrfunc(beta, nus1, nus2);
        roots = transpose(roots);
        roots = round(roots, 2);
        allroots(:, :, i) = transpose(roots);
        rootsun = unique(roots, 'rows');
        noroots(i, j) = length(rootsun);
    end
end
toc

%% Plot 3d
figure;
%stem3(betas, nus, noroots, 'LineStyle', 'none');
surf(betas, nus, noroots, 'edgeColor', 'none');
colormap turbo
xlabel('\beta')
ylabel('\nu')
zlabel('Roots')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Surface showing the number of roots')

%% Plot 2d
figure;
plot(betas, allroots(1,:,:));
colormap turbo
xlabel('\beta')
ylabel('Roots')
set(gca,'xscale','log')
title('Plot showing the number of roots')

%% Plot roots
scatter(roots(:, 1), roots(:, 2));