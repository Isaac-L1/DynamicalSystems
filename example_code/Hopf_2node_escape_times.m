% stochastic simulation of two conncted Hopf nodes - stochastic Heun method
% finding mean escape time for a variety of beta values
% close all
clear
 
%% set up variables

% model parameters
alpha = 0.1;           % noise amplitude
lambda = 0.8;           % excitability 
nu = 1-lambda;
omega = 0; 

bnum = 10;              % number of graph points to generate.
betavals = logspace(-3,2, bnum); % bnum log spaced points from 10^-3 and 10^2

% parameters for the numerical method
N = 10000000;           % number of steps to take
T = 10000;              % maximum time
h = 1e-3;                % time step
t = (0:h:T);            % t is the (discretised time) vector [0 1h 2h 3h ... Nh]

kmax=100;              % total number of simulations for each beta value to compute mean.

first_E = zeros(1,bnum);   % holder for first mean escape times
second_E = zeros(1,bnum);  % holder for second mean escape times
both_E = zeros(1,bnum);    % holder for both mean escape times

thresh= sqrt(1-sqrt(lambda))+0.1; % position of the limit cycle + a bit 

% set up timing
tic

%% computation

for j=1:bnum % find mean escape time for each beta value
    
    beta=betavals(j);
    fprintf(['\n Computing escape times for beta = ' num2str(beta) '\n']);
    
    Ef=zeros(1,kmax);         % set up escape time vector
    Es=zeros(1,kmax);         % set up escape time vector
    Eb=zeros(1,kmax);         % set up escape time vector
    
    for k=1:kmax
        y1 = zeros(size(t)); % prepare place to store realisation
        y2 = zeros(size(t));
        y1(1)=0;             % initial values
        y2(1)=0;
        t=0;
        n=1;
        
        while 1         % start taking steps
            
            % evaluate slope of deterministic bit left side of interval
            f1l=(lambda-1+1i*omega)*y1(n) + ...
                2*y1(n)*(abs(y1(n))^2) - y1(n)*(abs(y1(n))^4)+ ...
                beta*(y2(n) - y1(n));
            
            f2l=(lambda-1+1i*omega)*y2(n) + ...
                2*y2(n)*(abs(y2(n))^2) - y2(n)*(abs(y2(n))^4)+ ...
                beta*(y1(n) - y2(n));
            
            % prediction euler's step
            y1bar = y1(n) + h*f1l + alpha*sqrt(h)*(randn+sqrt(-1)*randn);
            y2bar = y2(n) + h*f2l + alpha*sqrt(h)*(randn+sqrt(-1)*randn);
            
            % correction
            f1r = (lambda-1+1i*omega)*y1bar + ...
                2*y1bar*(abs(y1bar)^2) - y1bar*(abs(y1bar)^4)+ ...
                beta*(y2bar - y1bar);
            
            f2r = (lambda-1+1i*omega)*y2bar + ...
                2*y2bar*(abs(y2bar)^2) - y2bar*(abs(y2bar)^4)+ ...
                beta*(y1bar - y2bar);
            
            % combine to find next step
            y1(n+1)= y1(n) + h*(f1l+f1r)/2 ... 
                + alpha*sqrt(h)*(randn+sqrt(-1)*randn);

            y2(n+1)= y2(n) + h*(f2l+f2r)/2 ...  
                + alpha*sqrt(h)*(randn+sqrt(-1)*randn);
            
            t=t+h;
            
            if (abs(y1(n+1))>thresh || abs(y2(n+1))>thresh) && Ef(k)==0  % first node has escaped
                Ef(k)=t;
            end
            
            if (abs(y1(n+1))>thresh && abs(y2(n+1))>thresh)  % both nodes have escaped
                Eb(k)=t;
                Es(k)=Eb(k)-Ef(k);        % escape time
                break
            end
            
            n=n+1;
        end


    end

    first_E(j) = mean(Ef);
    second_E(j) = mean(Es);
    both_E(j) = mean(Eb);

end

% stop timer
toc

%% plot it all
figure;hold on;

% plot the means
%plot(betavals,both_E,'k-','linewidth',2)
plot(betavals,first_E,'k--','linewidth',2) 
plot(betavals,second_E,'k:','linewidth',2)

legend('both','first', 'second','Location','SouthWest')
xlabel('\beta')
ylabel('E[\tau]','Rotation',0)
set(gca,'xscale','log');
box on

tnam=sprintf('%d beta values, %d simulations, lambda=%1.2f, alpha=%1.2f, h=%1.4f',bnum,kmax,lambda,alpha,h);
title(tnam);

