clear

%%
X = linspace(-0.6, 1.4, 301);
v = linspace(0, 0.1, 301);
V = zeros(length(v), length(X));

%%

for i = 1:length(v)
    for j = 1:length(X)
        vi = v(i);
        Xi = X(j);
        V(j, i) = (1/4)*Xi^4 - (1/3)*Xi^3 - (1/2)*Xi^2*vi + Xi*vi;
    end
end

%%
surf(v, X, V, 'edgeColor', 'none')
xlabel('\nu')
ylabel('X')
zlabel('V(X)')
title('Surface showing the Potential Landscape')

%%

X = linspace(-0.6, 1.4, 301);
nus = [0 0.001 0.01 0.05 0.1];
V = zeros(length(X), length(nus));

for i = 1:length(nus)
    for j = 1:length(X)
        vi = nus(i);
        Xi = X(j);
        V(j, i) = (1/4)*Xi^4 - (1/3)*Xi^3 - (1/2)*Xi^2*vi + Xi*vi;
    end
end

plot(X, V);
title('Potential Landscapes against \nu');
xlabel('x');
ylabel('V(x)');
hleg = legend('0', '0.001', '0.01', '0.05', '0.1','Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','\nu');
