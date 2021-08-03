clear


%% Variables

% model parameters 
alphas = [0.1, 0.125, 0.15]; % noise amplitudes
nu = logspace(-3, -1.5, 10);           % excitability 


% parameters for the numerical method
N = 10000000;           % number of steps to take
T = 10000;              % maximum time
h = 1e-3;                % time step
t = (0:h:T);            % t is the (discretised time) vector [0 1h 2h 3h ... Nh]

kmax=10000;              % total number of simulations for each nu value to compute mean.

E = zeros(3, 10);   % holder for mean escape times
S = zeros(3, 10);   %holder for stds


%% Simulation
tic

for b = 1:length(alphas)
    
    alpha = alphas(b);

    for i = 1:length(nu)
        thresh= 0.6; % position of the limit cycle + a bit 
        fprintf(['\n Computing escape times for nu = ' num2str(nu(i)) '\n']); %Print which step it is on

        En = zeros(1, kmax); %Set up escape time vector
        nui = nu(i);

        parfor j = 1:kmax %Find kmax realisations with escape times.
            x = zeros(size(t)); %Prepare vector to store realisations
            x(1) = -sqrt(nui); %Initial conditions
            tx = 0;
            n = 1;

            while 1 %Start taking time steps.
                %Evaluate deterministic part using Heun's scheme.
                k1 = F(x(n), nui);
                k2 = F((x(n)+h*k1 + alpha*randn*sqrt(h)), nui);
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

        E(b, i) = mean(En);
        S(b, i) = std(En);
    end



end

toc

%% Plot
figure;hold on;

% plot the means
plot(nu, E);

hleg = legend('0.1','0.125','0.15','Location','NW');
htitle = get(hleg,'Title');
set(htitle,'String','\alpha')
xlabel('\nu')
ylabel('E[\tau]','Rotation',0)
title('Mean escape times for a single node against \nu for different \alpha')
set(gca,'xscale','log');
box on
hold off

figure;hold on;

% plot the means
plot(nu, S);

hleg = legend('0.1','0.125','0.15','Location','NW');
htitle = get(hleg,'Title');
set(htitle,'String','\alpha')
xlabel('\nu')
ylabel('SD[\tau]','Rotation',0)
title('SD of escape times for a single node against \nu for different \alpha')
set(gca,'xscale','log');
box on


%% Function Declarations

function [x1] = F(x0, v)
    
    x1 = -(x0^2 - v)*(x0 - 1);

end