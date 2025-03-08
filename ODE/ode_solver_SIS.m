% Parameters
beta = 0.1;
alpha = 1/14;
%N = 245000;
N=10000;
cumI0 = 5;

% initial condition
X0 = [N-cumI0,cumI0]; % [S,I]

TT=0:0.1:60;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% initialization matriz to save data
XX=zeros(length(TT),length(X0));
XX(1,:) = X0;

for ii = 1:length(TT)-1
    % Call the ODE solver
    [t,px] = ode15s(@SIS, [TT(ii) TT(ii+1)], X0', options, beta, alpha, N);
    fprintf('time = %f \n',TT(ii+1));  
    X0 = px(end,:)';
    XX(ii+1,:) = X0;     
end

figure
hold on
plot(TT,XX(:,1),'b-','LineWidth',1.5)
plot(TT,XX(:,2),'r-','LineWidth',1.5)
xlabel('Days','Interpreter','latex')
ylabel('Persons','Interpreter','latex')
legend('S','I')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
hold off


figure
hold on
plot(TT,XX(:,2),'k-','LineWidth',1.5)
xlabel('Days','Interpreter','latex')
ylabel('Infected Persons','Interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
hold off