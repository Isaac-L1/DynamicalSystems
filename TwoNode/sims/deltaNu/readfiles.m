files = dir('*.csv');

for i=1:length(files)
    load(files(i).name, '-ascii');
end

%%


%%

Exn = simsdNs1.';


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

Data1 = [-0.001 "2->1" PY EEDY(1) SEDY(1) CEDY(1) EEDY(2) SEDY(2) CEDY(2);
         -0.001 "1->2" PX EEDX(1) SEDX(1) CEDX(1) EEDX(2) SEDX(2) CEDX(2)];
     
%%

dndata = [Data1;Data2;Data3;Data4;Data5;Data6];

dndata(:, 2) = [2;1;2;1;2;1;2;1;2;1;2;1];
for i = 1:12
    for j = 1:9
        dndata(i,j) = str2double(dndata(i,j));
    end
end
dndata = cast(dndata, 'double');

dnTable = array2table(dndata, 'VariableNames', {'deltaNu', 'FirstNode', 'P(X)', 'E(T1)', 'S(T1)', 'C(T1)', 'E(T2)', 'S(T2)', 'C(T2)'});

%%
data = table2cell(dnTable);
uitable('Data',data,'ColumnName',dnTable.Properties.VariableNames,...
    'RowName',dnTable.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


