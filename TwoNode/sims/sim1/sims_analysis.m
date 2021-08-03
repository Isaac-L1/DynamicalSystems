clear

load('escapetimes.mat');

nus = [0.001, 0.0025, 0.005, 0.01, 0.02, 0.03];
betas = [0, 0.002, 0.02];

Eminn = min(Exn, Eyn);
Emaxn = max(Exn, Eyn);
EDn = (Emaxn-Eminn);
indx = (Eminn == Eyn);
EminXn = Eminn;
EminXn(indx) = NaN;
EDXn = EDn;
EDXn(indx) = NaN;
indy = (Eminn == Exn);
EminYn = Eminn;
EminYn(indy) = NaN;
EDYn = EDn;
EDYn(indy) = NaN;

Ex = mean(Exn, 4);
Ey = mean(Eyn, 4);
Emin = mean(Eminn, 4);
Emax = mean(Emaxn, 4);
ED = mean(EDn, 4);
EDX = mean(EDXn, 4, 'omitnan');
EDY = mean(EDYn, 4, 'omitnan');
EminX = mean(EminXn, 4, 'omitnan');
EminY = mean(EminYn, 4, 'omitnan');

Sx = std(Exn, 0, 4);
Sy = std(Eyn, 0, 4);
Smin = std(Eminn, 0, 4);
Smax = std(Emaxn, 0, 4);
SD = std(EDn, 0, 4);
SDX = std(EDXn, 0, 4, 'omitnan');
SDY = std(EDYn, 0, 4, 'omitnan');
SminX = std(EminXn, 0, 4, 'omitnan');
SminY = std(EminYn, 0, 4, 'omitnan');

Cx = Sx./Ex;
Cy = Sy./Ey;
Cmin = Smin./Emin;
Cmax = Smax./Emax;
CD = SD./ED;
CDX = SDX./EDX;
CDY = SDY./EDY;
CminX = SminX./EminX;
CminY = SminY./EminY;

%%

for i = 1:3
    figure;
    pl = tiledlayout(2, 3, 'TileSpacing', 'tight');
    pln = 1;

    for j = 1:6
        nexttile(pln)
        hold on;
        plot(nus, reshape(Cmin(i,j,:), [1 6]));
        plot(nus, reshape(CD(i,j,:), [1 6]));
        %plot(nus1, Cmax(i,:));
        set(gca,'xscale','log');
        xlabel('\nu 1');
        ylabel('Coefficient of Variation');
        title(['\nu 2 = ' num2str(nus(j))]);
        %hleg = legend('First escape', 'Domino time', 'Location', 'NW');
        hold off;
        pln = pln + 1;
    end
    
    title(pl, ['Coefficients of variation for \beta = ' num2str(betas(i))]);

end

%%
figure; hold on;
histogram(reshape(EminYn(1, 4, 4, :), [1 10000]))