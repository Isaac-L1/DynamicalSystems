clear

%% Variables

% model parameters
alpha = 0.1;           % noise amplitude
nu = 0.005;           % excitability 


% parameters for the numerical method
N = 100000;           % number of steps to take
T = 100;              % maximum time
h = 1e-4;                % time step
t = (0:h:T);            % t is the (discretised time) vector [0 1h 2h 3h ... Nh]

kmax=4;              % total number of simulations for each nu value to compute mean.

E = zeros(1,10);   % holder for mean escape times

%% Simulation
tic

thresh= 0.6; % position of the limit cycle + a bit 
fprintf(['\n Computing escape times for nu = ' num2str(nu) '\n']); %Print which step it is on
    
En = zeros(1, kmax); %Set up escape time vector
th = zeros(1, length(t));
th(1) = thresh;

figure; hold on

for j = 1:kmax %Find kmax realisations with escape times.
    x = zeros(1, length(t)); %Prepare vector to store realisations
    
    x(1) = -sqrt(nu); %Initial conditions
    tx = 0;
    n = 1;
        
    while n < length(t) %Start taking time steps.
        %Evaluate deterministic part using Heun's scheme.
        k1 = F(x(n), nu);
        k2 = F((x(n)+h*k1 + alpha*randn*sqrt(h)), nu);
        xdet = x(n) + (h/2)*(k1+k2);
            
        %Include the stochastic value.
        x(n+1) = xdet + alpha*randn*sqrt(h);
        th(n+1) = thresh;
            
        tx = tx + h;
            
            
        %Step forward
        n = n + 1;
    end
    plot(t, x)
end

plot(t, th)
hold off

toc

%% Function Declarations

function [x1] = F(x0, v)
    
    x1 = -(x0^2 - v)*(x0 - 1);

end
