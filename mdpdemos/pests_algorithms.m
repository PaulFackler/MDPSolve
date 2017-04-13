 clc 
 disp('demonstrates three alternative algorithms for infinite horizon problems')
% parameters
 D = [0; 5; 20];  % damage costs
 C = 10;          % spraying cost
 delta = 0.95;    % discount factor

% state/action combinations
% columns: 1) value of action  2) value of state
 X = [1 1;
      1 2;
      1 3;
      2 1;
      2 2;
      2 3];
    
% rewards
 R = -[D D+C];
 
% probabilities with no action
 P1=[0.65  0.15 0.05;
     0.25  0.40 0.20;
     0.10  0.45 0.75];
% probabilities with spraying
 P2=[0.85 0.45 0.35;
     0.15 0.50 0.50;
     0    0.05 0.15];
 P=[P1 P2];

% set up model structure
 clear model
 model.name = 'Simple pest infestation example';
 model.X = X;
 model.R = R;
 model.d = delta;
 model.P = P;
 model.T = inf;
 options=struct('algorithm','p','print',1);
 results1 = mdpsolve(model,options);
 
 options=struct('algorithm','f','print',1);
 results2 = mdpsolve(model,options);
 
 options=struct('algorithm','f','print',1,'modpol',50);
 results3 = mdpsolve(model,options);
 
 disp('State and optimal Control using three algorithms')
 disp([results1.Xopt(:,[2 1]) results2.Xopt(:,[1]) results3.Xopt(:,[1])])
 disp('error in value function relative to policy iteration solution')
 fprintf('%1i  %10.6e  %10.6e\n',[(1:3)' results1.v-results2.v results1.v-results3.v]')

