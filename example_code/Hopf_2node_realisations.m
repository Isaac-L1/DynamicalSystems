% stochastic simulation of two conncted Hopf nodes - stochastic Heun method
% 2 nodes bidirectional coupling with varying strength
close all
clear

N = 2000000;             % number of steps to take
h = 1e-3;              % time step
T = N*h;                 % maximum time
t = (0:h:T);          % t is the (discretised time) vector [0 1h 2h 3h ... Nh]

% parameters for the model
alpha = 0.05;              % noise amplitude
lambda = 0.8;       % excitability
nu = 1-lambda;
omega = 0;               % frequency

BETVALS = [0.01, 0.1, 1];

thresh= sqrt(1-sqrt(lambda))+0.1; % position of the limit cycle + a bit 

for j=1:3
    beta = BETVALS(j);

    % plot it
    figure(j);hold on;
    
    % threshold
    plot([thresh thresh],[0 1.5],'k-', 'linewidth',1);
    plot([0 1.5],[thresh thresh],'k-', 'linewidth',1); 
    
    
    %% plot realisation
    
    y1 = zeros(length(t),1);   % prepare place to store locations
    y2 = zeros(length(t),1);   % prepare place to store locations
    for n=1:N           % start taking steps
        % evaluate slope of deterministic bit left side of interval
        f1l = (lambda-1+1i*omega)*y1(n) + 2*y1(n)*(abs(y1(n))^2) - y1(n)*(abs(y1(n))^4) + ...
                beta*(y2(n) - y1(n));
        f2l=(lambda-1+1i*omega)*y2(n) + 2*y2(n)*(abs(y2(n))^2) - y2(n)*(abs(y2(n))^4)+ ...
                beta*(y1(n) - y2(n));

        % prediction euler's step
        y1bar = y1(n) + h*f1l + alpha*sqrt(h)*(randn+sqrt(-1)*randn);
        y2bar = y2(n) + h*f2l + alpha*sqrt(h)*(randn+sqrt(-1)*randn);
        
        % correction
        f1r =  (lambda - 1 + 1i*omega)*y1bar + 2*y1bar*(abs(y1bar)^2) - y1bar*(abs(y1bar)^4) + ...
                beta*(y2bar - y1bar);
        f2r =  (lambda - 1 + 1i*omega)*y2bar + 2*y2bar*(abs(y2bar)^2) - y2bar*(abs(y2bar)^4) + ...
                beta*(y1bar - y2bar);
        
        % next y values
        y1(n+1) = y1(n) + h*(f1l+f1r)/2 + alpha*sqrt(h)*(randn+sqrt(-1)*randn);
        y2(n+1) = y2(n) + h*(f2l+f2r)/2 + alpha*sqrt(h)*(randn+sqrt(-1)*randn);

    end
  
    % find the radius
    r1 = abs(y1);
    r2 = abs(y2);
    
    % plot only some of the points or figure becomes very large
    k=10000;
    ind = 1:round(length(r1)/k,0):length(r1); % gives k points.
    plr1 = r1(ind);
    plr2 = r2(ind);

    
    if beta==0.01      
        plot(plr1,plr2,'-','color',[0 0.45 0.74],'markersize',10) % paper figure

        % equilibria
        pr = 0.081835174126658; %always there sink
        plot(pr,pr,'k.','markersize',30);

        pr = 1.376515855136776; %always there sink
        plot(pr,pr,'k.','markersize',30);
    
        pr = 0.3138584152651;  %undergoes pitchfork source
        plot(pr,pr,'square k','markersize',10,'markerfacecolor','k');

        pr1 = 1.374675430924813;
        pr2 = 0.132315831162051; %saddle node with sink
        plot(pr1,pr2,'k.','markersize',30);
        plot(pr2,pr1,'.k','markersize',30);

        pr1 = 1.374883076399233; %this equilibria saddle
        pr2 = 0.272090175109030;
        plot(pr1,pr2,'^k','markersize',10,'markerfacecolor','k');
        plot(pr2,pr1,'^k','markersize',10,'markerfacecolor','k');

        pr1 = 0.088893646175044; % two  more saddles
        pr2 = 0.320657958286937;
        plot(pr1,pr2,'^k','markersize',10,'markerfacecolor','k');
        plot(pr2,pr1,'^k','markersize',10,'markerfacecolor','k');
        
    elseif beta==0.1
        plot(plr1,plr2,'-','color',[0.47 0.67 0.19],'markersize',10) % paper figure

        % equilibria
        pr = 0.081835174126658; %always there sink
        plot(pr,pr,'k.','markersize',30);

        pr = 1.376515855136776; %always there sink
        plot(pr,pr,'k.','markersize',30);
    
        pr = 0.3138584152651;  %undergoes pitchfork source
        plot(pr,pr,'square k','markersize',10,'markerfacecolor','k');
    
        pr1 = 0.181294404650843;
        pr2 = 0.357718620184317;
        plot(pr1,pr2,'^k','markersize',10,'markerfacecolor','k');
        plot(pr2,pr1,'^k','markersize',10,'markerfacecolor','k');
        
    else
        plot(plr1,plr2,'-','color',[0.87 0.49 0],'markersize',10) % paper figure

        % equilibria
        pr = 0.081835174126658; %always there sink
        plot(pr,pr,'k.','markersize',30);

        pr = 1.376515855136776; %always there sink
        plot(pr,pr,'k.','markersize',30);
    
        pr = 0.3138584152651;  %undergoes pitchfork saddle
        plot(pr,pr,'^k','markersize',10,'markerfacecolor','k');
    end
    
    xlabel('r1')
    ylabel('r2','rotation',0)
    box on
    axis([0 1.5 0 1.5])
    
    set(gca,'YTick',0:0.5:1.5, 'XTick',0:0.5:1.5);
    
   

end


