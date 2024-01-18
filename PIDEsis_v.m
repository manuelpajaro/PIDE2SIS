%Parameters
cumI0 = 33; % initial number of infected persons
beta1 = 0.14; % infection rate 1
beta2 = 0.03; % infection rate 2
Tchange = 3; % day in which the infection rate changes

% Fixed parameters
N = 10000; 
b=1; 
gammaI=1/14; 

      
% Spatial discretization
xmin = 0;
xmax = 500;
deltax = 0.5;
x = xmin:deltax:xmax;

% Time definition
t0     = 0;
tmax   = 7;
Tsave  = tmax*10;
nt     = Tsave*10 + 1;
deltat = (tmax-t0)/(nt-1);
tl     = linspace(t0, tmax, nt);
fprintf('\n The time discretization is: %g \n',tl(2)-tl(1));

% Initial conditions (Gaussian density function)
PX0un = normpdf(x,cumI0,0.2);

% Normalized initial condition
PX0=PX0un/trapz(x,PX0un);

tic
% Computation of characteristics curves
xbar=x*exp(deltat*gammaI);
xbarlim=find(xbar >= x(end));
fprintf('\nThe length of X is: %g and there are %g points of Xbar biggest than Xmax \n',length(x),length(xbarlim));

% Initialization 
PX      = PX0;
Lix     = zeros(size(PX));

% Time Independent functions
e_x=exp(x/b);
e_lx=exp(-x/b);
e_x(isinf(e_x)==1)=realmax;
e_lx(isinf(e_lx)==1)=realmax;


% find the time where beta changes:
indTc_aux = find(tl<=Tchange);
indTc = indTc_aux(end);
% Input function, f(x), construction:
fx1 = beta1*x.*(N-x)/N; % times beta1 if t<indTc or beta2 if t>=indTc
fx2 = beta2*x.*(N-x)/N;
hx = 1-exp((x-N)/b);

% Saving only Tsave simulated times 
nt_sol = Tsave +1;
PX_sol      = zeros(nt_sol,length(x));
TT          = zeros(nt_sol,1);
TT(1)       = tl(1);    
kk_sol = 1;
jjkk=(nt-1)/(nt_sol-1)+1:(nt-1)/(nt_sol-1):nt;
PX_sol(1,:) = PX;   

for j = 2:nt
    % PX_bar construction using interpolation 
    PX_bar = interp1(x,PX,xbar);
    PX_bar(isnan(PX_bar)==1)=0;

    % fx 
    if j<indTc
        fx = fx1;
    else
        fx = fx2;
    end
       
    % Integral term computation by numerical integration (for x<1 Lix = 0)
    IL0 = length(find(x<1));
    Lix(IL0:end) = 1/b*e_lx(IL0:end).*cumtrapz(x(IL0:end),e_x(IL0:end).*fx(IL0:end).*PX(IL0:end));

    % Explicit method    
    PX = (PX_bar+deltat*Lix)./(1-deltat*(gammaI-hx.*fx)); 
    
    % Zero boundary condition
    PX(end) = 0;
    
    % Normalization: int_xmin^xmax(PX)dx=1
    PX=PX/trapz(x,PX);
            
    % Saving the solutio n fo the current time step
    if j==jjkk(kk_sol)
        fprintf('Time = %f \n',tl(j))
        PX_sol(kk_sol+1,:) = PX;
        TT(kk_sol+1)       = tl(j);
        kk_sol = kk_sol+1;
    end
end
toc

% Plotting the solution at final time
figure
hold on
plot(x,PX,'k-','LineWidth',1.5)
xlabel('Infected')
ylabel('Probability')
hold off

figure
mesh(x,TT,PX_sol)
xlabel('Infected')
ylabel('days')
zlabel('P(x,t)')
xlim([0 50])
ylim([0 7])

meanI = zeros(size(TT));
std   = zeros(size(TT)); 
for i=1:length(TT)
    meanI(i) = trapz(x,x.*PX_sol(i,:));
    std(i) = sqrt(trapz(x,x.^2.*PX_sol(i,:)));
end

figure
hold on
plot(TT,meanI,'k-','LineWidth',1.5)
plot(TT,meanI+std,'k--','LineWidth',1.25)
plot(TT,meanI-std,'k--','LineWidth',1.25)
hold off
xlabel('Days')
ylabel('Infected: $\mu_I \pm \sigma_I $','Interpreter','latex')

[X,Y] = meshgrid(x,TT);
figure
hold on
contour(Y',X',PX_sol',300)
plot(TT,meanI,'k-','LineWidth',1.5)
ylim([0,150])


[X,Y] = meshgrid(x,TT(2:end));
figure
hold on
contour(Y',X',PX_sol(2:end,:)',300)
plot(TT(2:end),meanI(2:end),'k-','LineWidth',1.5)
ylim([0,150])

%% to plot real data
Ireal = [33,40,40,40,39,37,36,35];% Gondomar 2020 10 25 beta1=0.14, beta2=0.03, Tchange=3
%Ireal = [27,33,39,42,47,50,49,52]; % Gondomar 2021 01 10 beta1=0.22, beta2=0.03, Tchange=5
[X,Y] = meshgrid(x,TT);
figure
hold on
contour(Y',X',PX_sol',300)
plot(TT,meanI,'k-','LineWidth',1.5)
plot([0,1,2,3,4,5,6,7],Ireal,'rs','MarkerSize',8,'LineWidth',3)
legend('P(I)_{PIDE}','I\_mean_{PIDE}','I\_Real')
ylim([0,100])
xlabel('$t$ (days)','Fontsize',22,'Interpreter','latex')
ylabel('$I$','Fontsize',22,'Interpreter','latex')
set(gca,'FontSize',12,'TickLabelInterpreter','latex')
hold off


