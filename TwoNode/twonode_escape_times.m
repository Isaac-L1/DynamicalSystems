parpool(16);

%Set up system specific constants and parameters of variation.
alpha = 0.1;
nus1 = logspace(-3, -1.5, 10);
nus2 = logspace(-3, -1.5, 10);

kmax = 5000; %Number of realisations.

T = 200;
h = 1e-3;

Ex = zeros(length(nus1), length(nus2), kmax);
Ey = zeros(length(nus1), length(nus2), kmax);

thresh = 0.8; %Threshold for escape.
lthresh = 0.5;
tic

for b = 1:length(nus2)
    
    nu2 = nus2(b);

    for i = 1:length(nus1)

        nu1 = nus1(i);

        for j = 1:kmax %Find kmax realisations with escape times.
            
            nums = randn(4, 8*T/h);
            
            x = zeros(1, T/h); %Prepare vector to store realisations
            y = zeros(1, T/h);

            x(1) = -0.1; %Initial conditions
            y(1) = -0.1;
            tx = 0;
            n = 1;
            xesc = 0; %Flag to check which nodes have escaped.
            yesc = 0;
            
            
            while 1 %Start taking time steps.
                
                if n>length(nums)
                    nums = randn(4, 8*T/h);
                end
                
                %Evaluate deterministic part using Heun's scheme.
                kx1 = -(x(n)^2 - nu1)*(x(n) - 1) + beta*(y(n)-x(n));
                ky1 = -(y(n)^2 - nu2)*(y(n) - 1) + beta*(x(n)-y(n));

                %Apply noise
                hkx1 = x(n) + h*kx1 + alpha*nums(1, n)*sqrt(h);
                hky1 = y(n) + h*ky1 + alpha*nums(2, n)*sqrt(h);
                
                kx2 = -(hkx1^2 - nu1)*(hkx1 - 1) + beta*(hky1-hkx1);
                ky2 = -(hky1^2 - nu2)*(hky1 - 1) + beta*(hkx1-hky1);

                %Apply heun scheme and noise.
                x(n+1) = x(n) + (h/2)*(kx1+kx2) + alpha*nums(3, n)*sqrt(h);
                y(n+1) = y(n) + (h/2)*(ky1+ky2) + alpha*nums(4, n)*sqrt(h);
                
                %Step forward
                n = n + 1;
                tx = tx + h;

                %Check if x has escaped.
                if xesc == 0 && x(n) > thresh
                    Ex(i, b, j) = tx;
                    xesc = 1;
                end

                %Check if y has escaped.
                if yesc == 0 && y(n) > thresh 
                    Ey(i, b, j) = tx;
                    yesc = 1;
                end
                
                %Check if x has escaped.
                if xesc == 1 && x(n) < lthresh
                    error('');
                end

                %Check if y has escaped.
                if yesc == 1 && y(n) < lthresh 
                    
                end

                %Check if both have escaped.
                if xesc == 1 && yesc == 1
                    break
                end
                
                
            end
            
        end

    end
end

%Write them to csv files.
writematrix(Ex, 'mean_escape_times_node1.csv');
writematrix(Ey, 'mean_escape_times_node2.csv');
writematrix(Sx, 'sd_escape_times_node1.csv');
writematrix(Sy, 'sd_escape_times_node2.csv');
toc
