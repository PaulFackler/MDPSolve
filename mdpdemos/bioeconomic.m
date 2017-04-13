% Bioeconomic Model

% Enter model parameters
  T     = 10;                   % foraging periods
  emax  =  8;                   % energy capacity
  e     = [2  4  4];            % energy from foraging
  p     = [1.0 0.7 0.8];        % predation survival probabilities
  q     = [0.5 0.8 0.7];        % foraging success probabilities

% Construct state space
  S = (0:emax)';                % energy levels
  n = length(S);                % number of states
  m = 3;                        % number of actions

% Construct reward function
  R = zeros(n,m);

% Construct state transition probability matrix
  P = zeros(n,n,m);
  for k=1:m
    P(1,1,k) = 1;
    for j=2:n
       % does not survive predation
       snext = 0;           i=match(snext,S); P(i,j,k) = P(i,j,k) + (1-p(k)); 
       % survives predation, finds food
       snext = S(j)-1+e(k); i=match(snext,S); P(i,j,k) = P(i,j,k) + p(k)*q(k);
       % survives predation, finds no food
       snext = S(j)-1;      i=match(snext,S); P(i,j,k) = P(i,j,k) + p(k)*(1-q(k));
    end
  end
  P=reshape(P,n,n*m);
   
% Initialize terminal value function
  vterm = ones(n,1);            % terminal value: survive
  vterm(1) = 0;                 % terminal value: death
  
% Pack model structure
  clear model
  model.R      = R;
  model.P      = P;
  model.T      = T;
  model.d      = 1;
  model.vterm  = vterm;

% Solve finite-horizon model using backward recursion
  options=struct('keepall',1);
  results = mdpsolve(model,options);
  v=results.v;
  
% Plot survial probabilities, period 1
  figure(1); 
  h=bar(S,v(:,1),1);  set(h,'FaceColor',[.75 .75 .75])

  axis([-.5 emax+.5 0 1]);
  title('Survival Probability (Period 0)');
  xlabel('Stock of Energy'); ylabel('Probability');
  
% Plot survial probabilities, period 5
  figure(2); 
  h=bar(S,v(:,6),1);   set(h,'FaceColor',[.75 .75 .75])
  axis([-.5 emax+.5 0 1]);
  title('Survival Probability (Period 5)');
  xlabel('Stock of Energy'); ylabel('Probability');
  