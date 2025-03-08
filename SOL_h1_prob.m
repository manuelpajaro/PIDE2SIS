% Fixed parameters
N = 10000; 
b=1; 
gammaI=1/14; 

% beta
BETA = 0.001:0.001:0.1;
Pext = zeros(size(BETA));


for i=1:length(BETA)
    beta=BETA(i);
    
    % Number of infected
    x1=0.0001:0.0001:1;
    x2=2:1:5000;
    x=[x1,x2];
    
    % Solution (without normalization)
    Px_un = 1./x.*exp(-x/b+x/gammaI*beta-x.^2/(2*gammaI*N)*beta);
    Px = Px_un/trapz(x,Px_un);  
    
    %% Prob <=1
    Pext(i) = trapz(x1,Px(1:length(x1)));
    %Pper = trapz(x(length(x1):end),Px(length(x1):end))
    %Pext+Pper
end

%% R0 = (b*beta)/gammaI ==1
beta_star=gammaI/b;
sizeL = 14; % size for figures axis,labels,titles

figure
hold on
plot(BETA,Pext,'k-','LineWidth',1.5)
plot(beta_star*ones(1,11),0:0.1:1,'k:','LineWidth',1.5)
text(beta_star,0.85,'\leftarrow \beta^*','Fontsize',sizeL)
xlabel('$\beta$','Fontsize',sizeL,'Interpreter','latex')
ylabel('$P(I\leq 1)$','Fontsize',sizeL,'Interpreter','latex')
title(['b = ',num2str(b)],'Interpreter','latex')
set(gca,'FontSize',sizeL,'TickLabelInterpreter','latex')
hold off