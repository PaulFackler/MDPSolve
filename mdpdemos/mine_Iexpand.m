  clc
% Model parameters
  price = 1;                    % price of ore
  sbar  = 100;                  % initial ore stock
  delta = 0.9;                  % discount factor  

% Construct state and action spaces
  S = (0:sbar)';                % vector of states
  A = (0:sbar)';                % vector of actions
  X = rectgrid(S,A);            % matrix of state/action combinations
  X = X(X(:,2)<=X(:,1),:);      % eliminate infeasible actions
  
  Iexpand = X(:,1)-X(:,2)+1;
  Ix      = X(:,1)+1;
  
  P = speye(length(S));
  R = (price-X(:,2)./(1+X(:,1))).*X(:,2);
  
  clear model
  model.R  = R;
  model.P  = P;
  model.d  = delta;
  model.Iexpand = Iexpand;
  model.Ix = Ix;
  model.colstoch=1;
  options=struct('algorithm','p','getAopt',1);
  results = mdpsolve(model,options);
  if ~isempty(results.errors)
    mdpreport(results);
    return
  end
  v=results.v;
  Aopt=getA(X,results.Ixopt,2);
  
% Plot optimal value function
  figure(1); plot(S,v);
  title('Optimal Value Function');
  xlabel('Stock'); ylabel('Value');

% Plot optimal policy function
  figure(2); plot(S,Aopt); 
  title('Optimal Extraction Policy');
  xlabel('Stock'); ylabel('Extraction');
  
  %textable([X Iexpand Ix],[0 0 0 0],{'$S$','$A$','Iexpand','Ix'})
