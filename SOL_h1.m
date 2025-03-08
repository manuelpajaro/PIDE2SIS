% Fixed parameters
N = 1000; 
b=1; 
gammaI=1/14; 

% beta
%beta = 0.1; % infection rate
beta_star=gammaI/b;
beta=0.1;

% Number of infected
%x0=0.01:1:5000;
x1=0.0001:0.0001:1;
x2=2:1:5000;
xx=[x1,x2];

%%
x=xx;

% Solution 
Px_un = 1./x.*exp(-x/b+x/gammaI*beta-x.^2/(2*gammaI*N)*beta);
Px = Px_un/trapz(x,Px_un);

% Solution 
b1=2;
beta1=beta/b1;
Px_un1 = 1./x.*exp(-x/b1+x/gammaI*beta1-x.^2/(2*gammaI*N)*beta1);
Px1 = Px_un1/trapz(x,Px_un1);

% Solution 
b2=3;
beta2=beta/b2;
Px_un2 = 1./x.*exp(-x/b2+x/gammaI*beta2-x.^2/(2*gammaI*N)*beta2);
Px2 = Px_un2/trapz(x,Px_un2);

figure
plot(x,Px_un)

figure
plot(x,Px,'k-','LineWidth',1.5)
xlabel('I','Interpreter','latex')
ylabel('P(I)','Interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
%ylim([0,0.0001])
%xlim([0,1000])

%% figure
sizeL = 14; % size for figures axis,labels,titles
Istar = max(0,(1-gammaI/(b*beta))*N);

figure
hold on
plot(x,Px,'k-','LineWidth',1.5)
plot(x,Px1,'b-','LineWidth',1.5)
plot(x,Px2,'r-','LineWidth',1.5)
plot(Istar*ones(1,11),linspace(0,1.5*max(Px),11),'k:','LineWidth',1.5)
text(Istar,1.25*max(Px),'\leftarrow I^*','Fontsize',sizeL)
xlabel('$x$','Fontsize',sizeL,'Interpreter','latex')
ylabel('$P_s(x)$','Fontsize',sizeL,'Interpreter','latex')
legend('$b=1$','$b=2$','$b=3$','Interpreter','latex')
%title(['b = ',num2str(b)],'Interpreter','latex')
set(gca,'FontSize',sizeL,'TickLabelInterpreter','latex')
ylim([0,1.5*max(Px)])
xlim([0,1000])
hold off
