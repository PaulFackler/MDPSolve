
 disp('A simple example - managing a pest')

% parameters
 D = [0; 5; 20];  % damage costs
 C = 10;          % spraying cost
 delta = 0.95;    % discount factor
 
 reward=@(A,S) -D(S) - C*(A==1);

% probabilities with no action
 P1=[0.65  0.15 0.05;
     0.25  0.40 0.20;
     0.10  0.45 0.75];
% probabilities with spraying
 P2=[0.85 0.45 0.35;
     0.15 0.50 0.50;
     0    0.05 0.15];
 
 S=[1;2;3];
 A=[0;1];
 X=rectgrid(A,S);
 Stran=rvdef('d',[P1 P2],S);
 %D=add2diagram(D,name,type,obs,parents,cpd,loc,attachments)
 D=add2diagram([],'pest level', 's',1,{},S,[.3 .6]);
 D=add2diagram(D, 'treat',      'a',1, {},A,[.3 .4]);
 D=add2diagram(D, 'pest level+','f',1,{'treat','pest level'},Stran,[.8 .6]);
 D=add2diagram(D, 'utility',    'r',1,{'treat','pest level'},reward,[.6 .4]);
 
 figure(1); clf
 drawdiagram(D)
 
 model=d2model(D);
 model.d = delta;
 options=struct('algorithm','i');
 results = mdpsolve(model,options);
 
 disp('Simple Pest Control Problem')
 disp('State, Optimal Control and Value')
 fprintf('%1i  %1i  %8.2f\n',[results.Xopt(:,[2 1]) results.v]')
%%
T=25;
reps=50;
s0=2+zeros(reps,1);
keepall=1;
Aopt=X(results.Ixopt,1);
[Y,z]=dsim(D,s0,T,Aopt,[],[],keepall);