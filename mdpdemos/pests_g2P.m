% pest management problem
% uses g2P
 clc 
 disp('Demonstrates the use of g2P for a discrete problem')
 % create matrix of state/action combinations
 X=[ ...
   1 1
   2 1
   3 1
   1 2
   2 2
   3 2];
 D=[0;5;20];  % damage costs
 C=10;        % spraying cost
 delta=0.95;  % discount factor
 % probabilities with no action
 P1=[0.65  0.15 0.05;
     0.25  0.40 0.20;
     0.10  0.45 0.75];
 % probabilities with spraying
 P2=[0.85 0.45 0.35;
     0.15 0.50 0.50;
     0    0.05 0.15];
 P=[P1 P2];
 R=-[D D+C]; 
 [Ix,S]=getI(X,1);
 
 g = @(X,e) e+zeros(size(X,1),1);
 Pf=g2P(g,S,X,S,P);
 
 % set up model structure
 clear model
 model.X        = X;
 model.discount = delta;
 model.R        = R;
 model.P        = Pf;
 % call solution procedure
 options=struct('algorithm','p','getAopt',1);
 results=mdpsolve(model,options);
    
 disp('Simple Pest Control Problem')
 disp('State, Optimal Control and Value')
 fprintf('%2i  %2i %8.4f\n',[S results.Xopt(:,2) results.v]')
 
 