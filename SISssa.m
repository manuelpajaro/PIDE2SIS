% We must define the number of reactions (j), the number of species (n) the
% values of the j rate constants (c0) , the initial state (n0), and the 
% state change vector (nu) 
% times at which we return the number of species obtained by SSA moreover
% Tgrig(1) is the initial time considered to start SSA and Tgrig(end) is
% the last time or stop condition to SSA.
Tgrid=0:1:7; % days

reaction_number=2;
species_number=2;

N = 10000;
cumI0 = 40;

x0=[N-cumI0,cumI0]; %[SI]

% reaction asociated (2 reactions, 2 species)
% S --> Susceptible (S=N-I)
% I --> Infectious (species 1)
% N total population susceptible
%%%%%%%% 
% S --> I   : rc = beta*S*I/N
% I --> R   : rc = gammaI*I
        
beta = 0.1;
gammaI = 1/14;
        
c0=zeros(1,reaction_number);
c0=c0+[beta/N  gammaI];
nu0=zeros(reaction_number,species_number);
nu=nu0+[-1  1  ;
         1 -1 ]; % state change vectors dim jxn
propensity=@(x) [c0(1)*x(1)*x(2) c0(2)*x(2)];

nsimula=1e4;

tic
simulation=SSA_mpd(propensity,nu,x0,Tgrid,nsimula);
toc

% plot the nsim simulation
nsim=1;
for i=1:length(x0)
    figure
    hold on
    plot(Tgrid,simulation{nsim}(i,:),'k-','LineWidth',1.5)
    xlabel('Time','Fontsize',14,'Fontname','Times new Roman');
    ylabel(['Species ',num2str(i)],'Fontsize',14,'Fontname','Times new Roman');
    set(gca,'FontSize',14);
    hold off
end

if nsimula>1
    % histogram
    indthist=length(Tgrid); % time to plot
    for i=1:length(x0)
        % histogram in the last time, we can select other
        %indthist=length(Tgrid);
        thist=Tgrid(indthist);
        YThist=zeros(1,nsimula);
        for j=1:nsimula
            YThist(j)=simulation{j}(i,indthist);
        end
        xmax=max(YThist);
        xmin=min(YThist);
        barhist=20;%max(floor(xmax/10),40);
        xx=linspace(xmin,xmax,barhist);
        figure
        hist(YThist,xx)
        hh=hist(YThist,xx);

        xplot=hh/trapz(xx,hh);
        figure
        plot(xx,xplot,'k--','LineWidth',1.75)
        hold on
        xlabel(['Species ',num2str(i)],'Fontsize',14,'Fontname','Times new Roman');
        ylabel('Probability','Fontsize',14,'Fontname','Times new Roman');
        set(gca,'FontSize',14);
        hold off 
    end
end



