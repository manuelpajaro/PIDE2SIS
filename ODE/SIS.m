function dy = SIS(t,y,beta,alpha,N)
  S = y(1);
  I = y(2);
  dS = -beta/N*S*I+alpha*I;
  dI = beta/N*S*I - alpha*I;
  dy = [dS,dI]';
end