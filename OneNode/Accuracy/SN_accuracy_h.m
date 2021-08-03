clear

%% Variables

% model parameters 
alphas = 0.1; % noise amplitudes
nu = 0.01;           % excitability 


% parameters for the numerical method
T = 10000;              % maximum time
hn = logspace(0, -4, 20);                % time step           % t is the (discretised time) vector [0 1h 2h 3h ... Nh]

times = 20;
kmax=1000;              % total number of simulations for each nu value to compute mean.

E = zeros(1, 20);   % holder for mean escape times
S = zeros(1, 20);
Er = zeros(1, 20);
Er(1) = 0;



%% Simulation
tic

for tn = 1:length(hn)
    
    fprintf(['\n Computing escape times for h = ' num2str(hn(tn)) '\n']);
    
    h = hn(tn); %h is the step size
    t = (1:h:T); %t is the discretised time vector
    Et = zeros(1, times);

    for ti = 1:times

        for b = 1:length(alphas)

            alpha = alphas(b);

            for i = 1:length(nu)
                thresh= 0.6; % position of the limit cycle + a bit 
                %fprintf(['\n Computing escape times for nu = ' num2str(nu(i)) '\n']); %Print which step it is on

                En = zeros(1, kmax); %Set up escape time vector

                for j = 1:kmax %Find kmax realisations with escape times.
                    x = zeros(size(t)); %Prepare vector to store realisations
                    x(1) = -sqrt(nu(i)); %Initial conditions
                    tx = 0;
                    n = 1;

                    while 1 %Start taking time steps.
                        %Evaluate deterministic part using Heun's scheme.
                        k1 = F(x(n), nu(i));
                        k2 = F((x(n)+h*k1 + alpha*randn*sqrt(h)), nu(i));
                        xdet = x(n) + (h/2)*(k1+k2);

                        %Include the stochastic value.
                        x(n+1) = xdet + alpha*randn*sqrt(h);

                        tx = tx + h;

                        %Check if it has escaped.
                        if x(n+1) > thresh
                            En(j) = tx;
                            break
                        end

                        %Step forward
                        n = n + 1;
                    end
                end

                Et(b, ti) = mean(En);

            end
        end
    end

    E(b, tn) = mean(Et);
    S(b, tn) = std(Et);

end

toc

%% Plot
figure;hold on;

% plot the means
plot(hn, E);

%hleg = legend('0.1','0.125','0.15','0.175','0.2','Location','NW');
%htitle = get(hleg,'Title');
%set(htitle,'String','\alpha')
xlabel('h')
ylabel('E[\tau]','Rotation',0)
title('Escape times for a single node for different step sizes')
set(gca,'xscale','log');
box on
hold off

figure; hold on;
xlabel('h')
ylabel('SD[\tau]','Rotation',0)
title('SD of escape times for a single node for different step sizes')
plot(hn, S);
set(gca,'xscale','log');
box on
hold off

figure; hold on;
av = mean(E);
Ea = abs(E-av);
xlabel('h')
ylabel('Deviation','Rotation',0)
title('Deviation from mean of escape times for a single node for different step sizes')
plot(hn, Ea);
set(gca,'xscale','log');
box on
hold off

%% Function Declarations

function [x1] = F(x0, v)
    
    x1 = -(x0^2 - v)*(x0 - 1);

end