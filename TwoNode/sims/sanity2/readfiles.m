files = dir('*.csv');

for i=1:length(files)
    load(files(i).name, '-ascii');
end

%%

Exn = zeros(2, 10000);

Exn(:, 1:2000) = sims1;
Exn(:, 2001:4000) = sims2;
Exn(:, 4001:6000) = sims3;
Exn(:, 6001:8000) = sims4;
Exn(:, 8001:10000) = sims5;


%%
Exn = Exn.';

X1 = Exn(1,:);
X2 = Exn(2,:);
%X3 = Exn(3,:);
Emin = min(Exn);
Emax = max(Exn);

Ex1 = mean(Exn(1, :));
Ex2 = mean(Exn(2, :));
%Ex3 = mean(Exn(3, :));

Sx1 = std(Exn(1, :));
Sx2 = std(Exn(2, :));
%Sx3 = std(Exn(3, :));

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
%ED(3, :) = EO(3, :) - EO(2, :);

indX = (pos(1, :) == 1);
EDX = ED(:, indX);
indY = (pos(2, :) == 1);
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

Data = [EEDX SEDX CEDX EEDY SEDY CEDY];