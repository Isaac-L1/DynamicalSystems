files = dir('*.csv');

for i=1:length(files)
    load(files(i).name, '-ascii');
end

Exn = zeros(4, 10000);

Exn(:, 1:2000) = sims1;
Exn(:, 2001:4000) = sims2;
Exn(:, 4001:6000) = sims3;
Exn(:, 6001:8000) = sims4;
Exn(:, 8001:10000) = sims5;


%%

X1 = Exn(1,:);
X2 = Exn(2,:);
X3 = Exn(3,:);
X4 = Exn(4,:);
Emin = min(Exn);
Emax = max(Exn);

Ex1 = mean(Exn(1, :));
Ex2 = mean(Exn(2, :));
Ex3 = mean(Exn(3, :));
Ex4 = mean(Exn(4, :));

Sx1 = std(Exn(1, :));
Sx2 = std(Exn(2, :));
Sx3 = std(Exn(3, :));
Sx4 = std(Exn(4, :));

pos = zeros(4, 10000);
type = zeros(1, 10000);
ED = zeros(4, 10000);
EO = zeros(4, 10000);

% for i = 1:10000
%     mn = find(Exn(:, i) == Emin(i));
%     pos(mn, i) = 1;
%     mx = find(Exn(:, i) == Emax(i));
%     pos(mx, i) = 4;
%     
%     
%     
% %     if all(pos(:, i) == [3;2;1])
% %         type(i) = 1;
% %     elseif all(pos(:, i) == [3;1;2])
% %         type(i) = 2;
% %     elseif all(pos(:, i) == [2;3;1])
% %         type(i) = 3;
% %     elseif all(pos(:, i) == [2;1;3])
% %         type(i) = 4;
% %     elseif all(pos(:, i) == [1;2;3])
% %         type(i) = 5;
% %     elseif all(pos(:, i) == [1;3;2])
% %         type(i) = 6;
% %     end
%     
%     EO(pos(1, i), i) = Exn(1, i);
%     EO(pos(2, i), i) = Exn(2, i);
%     EO(pos(3, i), i) = Exn(3, i);
%     EO(pos(4, i), i) = Exn(4, i);
%     
% end

EO = sort(Exn);

ED(1, :) = EO(1, :);
ED(2, :) = EO(2, :) - EO(1, :);
ED(3, :) = EO(3, :) - EO(2, :);
ED(4, :) = EO(4, :) - EO(3, :);

EED = mean(ED, 2);
SED = std(ED, 0, 2);

CED = SED./EED;


