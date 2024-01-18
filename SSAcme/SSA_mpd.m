% X(t)=(X_1(t),...,X_N(t)) state vector 
% nu_j=(nu(1j),...,nu(Nj))state-change vector
% c_jdt probability that any particular S1 molecule will so react in the
% next infenitesimal time dt. (x1.c_jdt if we have x1  particles of S1)
% a_j(x)=c_jx1 propensity function
% a_j(x)=c_jx1x2 if x1 of S1 and x2 of S2
% a_j(x)=c_j1/2x1(x1-1) if S1+S1 -> Product
% Chemical Master Ecuation
% dP(x,t|x0,t0)/dt=sum_j=1^M[aj(x-nu_j)P(x-nu_j,t|x0,t0)-a_j(x)P(x,t|x0,t0)]

% p(tau,j|x,t)=a_j(x)exp(-a0(x)tau) where ao(x)=sum_j´=(x)^M[a_j´(x)]
% tau is an exponencial random variable with mean (and standard desviation)
% 1/a0(x) and j statistically independent integer ramdon variable with
% point probabilities aj(x)/a0(x)

% Let r1 and r2 random numbers from the uniform distribution in the unit
% interval, them:
% tau=1/a0(x)ln(1/r1)
% jj = the smallest integer satisfying sum_j=1^jj[aj(x)]>r2a0(x)

function simulation=SSA_mpd(propensity,nu,x0,Tgrid,nsimula)
    % Initialization
    t0=Tgrid(1); % inicial time
    Tend=Tgrid(end);
    maxlength = 3e5; % aprox number of reactions fire in [t0,Tend]
    simulation=cell(nsimula,1);

    for nsim=1:nsimula
    x=x0; 
    t=t0;  
    ng=1;

    YY=zeros(length(x0),maxlength);
    TTT=zeros(1,maxlength);
    TTT(1)=t;
    YY(:,1)=x0';
    while t< Tend,
        propensities=propensity(x);

        a0=sum(propensities);

        % generate r1 and r2 (uniform distribution in the unit interval)
        r1=random('unif',0,1);
        r2=random('unif',0,1);

        % values for tau and jj
        tau=1/a0*log(1/r1);

        jj=0;
        aa=0;
        for i=1:size(nu,1)
            aa=aa+propensities(i);
            if aa>r2*a0
                jj=i;
                break;
            end        
        end

        % next reaction
        if jj==0
            t=Tend;
%            fprintf('x = [%f %f %f %f %f %f %f] \n',x(1),x(2),x(3),x(4),x(5),x(6),x(7));
        else
            t=t+tau;
            x=x+nu(jj,:);
        end
       ng=ng+1;
    TTT(ng)=t; 
    YY(:,ng)=x';
    end 

    if TTT(end)==0
        maxt=max(TTT);
        indt=find(TTT==maxt);
        fprintf('Simulation number: %f, %g reactions fired \n',nsim,indt);
        TTT=TTT(1:indt);
        YY=YY(:,1:indt);
    else
        fprintf('Simulation number: %f, %g reactions fired \n',nsim,length(TTT));
    end


    XX=zeros(length(x0),length(Tgrid));
    for i=1:length(x0)
        XX(i,:)=interp1(TTT,YY(i,:),Tgrid,'nearest');
    end
    simulation{nsim}=XX;

    end

end










