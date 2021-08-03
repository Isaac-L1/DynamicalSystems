alpha = 0.1;

D = 0.2;
H = 0.0008;
kmax = 1000;
escapes = zeros(1, kmax);
h = 1e-3;
thresh = 0.8;

nu = (D^2)/4;
A = (6*H)/(D^3);

for j = 1:kmax %Find kmax realisations with escape times.
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
            escapes(j) = tx;
            break
        end

        %Step forward
        n = n + 1;
        x0 = x1;
    end
end

meanEscapes = mean(escapes);
