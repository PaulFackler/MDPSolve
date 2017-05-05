  clc
  disp('Mine management')
  delta = 0.95;                   % discount rate
  p     = 1;                      % ore price
  c     = 0.5;                    % extraction cost parameter
  Sbar  = 100;                    % maximum ore level
  C = @(s,a) c*a.^2./(1e-300+s);  % extraction cost function
  g = @(s,a) s-a;                 % state transition function
  
  n     = 101;                    % number of ore levels
  sigma = Sbar/(n-1);             % extraction increment
  s     = sigma*(0:(n-1))';       % state values
  a     = s;                      % action values
  
  X=zeros(n*(n+1)/2,2);           % allocate memory for state/action combinations
  Ix=zeros(n*(n+1)/2,1);          % vector of state indices
  k=0;
  for i=1:n
    for j=1:n
      if a(j)<=s(i)
        k=k+1;
        X(k,:)=[s(i) a(j)];
        Ix(k)=i;
      end
    end
  end
  S=X(:,1);
  A=X(:,2);
  P = speye(n);

  R=p*A-C(S,A);
  clear model
  model.P = P;
  model.R = R;
  model.d=delta;
  model.Ix=Ix;
  model.Iexpand=match(S-A,s);
  model.colstoch=true;
  results = mdpsolve(model);
  v=results.v; x=results.Ixopt;
  a=A(x);
  
  alpha=roots([delta^2 -(2*delta*p+4*(1-delta)*c) p^2]);
  alpha=min(alpha);
  beta=(p-delta*alpha)/(2*c);
  
   
  figure(1)
  plot(s,v,'k')
  xlabel('S')
  ylabel('V(S)')
  
  figure(2)
  plot(s,[a,beta*s],'k')
  xlabel('S')
  ylabel('A^*(S)')
  
  figure(3)
  plot(s,v-alpha*s,'k')
  xlabel('S')
  ylabel('value function error')
  
  disp('approximate and exact alpha')
  disp([s\v alpha])
  disp('approximate and exact beta')
  disp([s\a beta])