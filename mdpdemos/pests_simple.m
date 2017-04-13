 clc 
 disp('A simple example - managing a pest')

% parameters
 D = [0; 5; 20];  % damage costs
 C = 10;          % spraying cost
 delta = 0.95;    % discount factor

% state/action combinations 
% columns: 1) value of action (no spray/spray) 
%          2) value of state (infestation level)
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
 results = mdpsolve(model);
 
 disp('Simple Pest Control Problem')
 disp('State, Optimal Control and Value')
 fprintf('%1i  %1i  %8.2f\n',[results.Xopt(:,[2 1]) results.v]')


