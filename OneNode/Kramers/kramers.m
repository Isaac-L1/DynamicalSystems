alpha = 0.1;

distances = linspace(0.02, 0.2, 20);
heights = linspace(1e-3, 2e-2, 40);
kmax = 5000;
escapes = zeros(length(distances), length(heights), kmax);
meanEscapes = zeros(length(distances), length(heights));
h = 1e-3;

%%
tic
for di = 1:length(distances)
    
    D = distances(di);
    
    for hi = 1:length(heights)
        
        H = heights(hi);
        
        nu = (D^2)/4;
        A = (6*H)/(D^3);
        thresh = sqrt(nu);
        
        parfor j = 1:kmax %Find kmax realisations with escape times.
            x0 = -sqrt(nu); %Initial conditions
            tx = 0;
            n = 1;

            while 1 %Start taking time steps.
                %Evaluate deterministic part using Heun's scheme.
                k1 = -(x0^2 - nu)*(A*x0 - A);
                
                hk1 = h*k1 + alpha*randn*sqrt(h);
                
                k2 = -((x0+hk1)^2 - nu)*(A*(x0+hk1) - A);
                
                x1 = x0 + (h/2)*(k1+k2) + alpha*randn*sqrt(h);

                tx = tx + h;

                %Check if it has escaped.
                if x1 > thresh
                    escapes(di, hi, j) = tx;
                    break
                end

                %Step forward
                n = n + 1;
                x0 = x1;
            end
        end
        
        meanEscapes(di,hi) = mean(escapes(di, hi, :));
        
    end
end
toc

%%

expEsc = zeros(length(distances), length(heights));
for i = 1:length(distances)
    for j = 1:length(heights)
        D = distances(i);
        H = heights(j);
        nu = D^2 / 4;
        A = (6*H) / D^3;
        expEsc(i, j) = pi / (A*sqrt(nu*(1-nu))) * exp((8*A*nu^(3/2)) / (3 * alpha^2));
    end
end

%%

exps = expEsc(:,21:40);
means = meanEscapes(:,21:40);
pref = sum(sum(means))/sum(sum(exps));

%%

figure(4); hold on;
surf(heights, distances, meanEscapes, 'edgecolor', 'none');
surf(heights, distances, expEsc, 'edgecolor', 'none');
xlabel('Heights');
ylabel('Distances');
zlabel('E( \tau )');
% set(gca, 'xscale', 'log');
% set(gca, 'yscale', 'log');
% set(gca, 'zscale', 'log');
view(3);

%%

figure(4); hold on;
plot(heights, meanEscapes(20,:));
plot(heights, expEsc(20,:));
hleg = legend('Computed','Approximation','Location','NW');
title('Mean escapes compared with approximation');
xlabel('Height');
ylabel('E(\tau)');

%% Chi-Square

sumCols = sum(meanEscapes);
sumRows = sum(meanEscapes, 2);
total = sum(sumCols);
expMeans = sumRows*sumCols / total;
chisqr = (meanEscapes - expMeans).^2;
chis = sum(sum(chisqr));
degs = (length(distances)-1)*(length(heights) - 1);

P = gammainc(degs/2, chis/2)/gamma(degs/2);

%Critical value at 10% significance is 790.744
%There is not enough evidence to reject the null hypothesis.





