% pests_pomdp Pest infestation control with monitoring
% Based on HaightPolasky10 
close all
clc
p=200;        % number of belief states
Xvals=[ ...
   1 1
   2 1
   3 1
   1 2
   2 2
   3 2];
delta=0.95;   % discount factor
D=[0;10;20];  % damage costs
Cm=inf;         % monitoring cost
Ct=20;        % treatment cost
T=20;         % time horizon
 
Pn=[0.8 0.0 0.0
    0.2 0.8 0.0
    0.0 0.2 1.0];
 
Pt=[0.9 0.8 0.6
    0.1 0.2 0.4
    0.0 0.0 0.0];
  
Qn=[0.5 0.3 0.1
    0.5 0.4 0.4
    0.0 0.3 0.5];
  
Qm=[1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0];
 
P=[Pn Pt];
R=-[D D+Ct];

% solve perfect certainty problem
model=struct('P',P,'R',R,'discount',delta,'T',T);
options=struct('algorithm','p','maxit',600);
results=mdpsolve(model,options);
vc=results.v; xc=results.Ixopt; 

% solve the POMDP
P=[Pn Pn Pt Pt];
Q=[Qn Qm Qn Qm];
R=-[D D+Cm D+Ct D+Cm+Ct];
options=struct('Qtype',1,'Rtype',2);
[b,Pb,Rb]=pomdp(p,P,Q,R,options);

model=struct('P',Pb,'R',Rb,'discount',delta,'T',T);
options=struct('algorithm','p','maxit',600,'nochangelim',500,'prtiters',0);
results=mdpsolve(model,options);
v=results.v; x=results.Ixopt; 

% display certainty results
Avals=repmat((1:4),size(b,1),1); 
ii=[find(b(:,1)==1);find(b(:,2)==1);find(b(:,3)==1)];
disp('State, certainty value and action, POMDP value and action (when state is known with certainty)')
fprintf('%1i  %10.4f  %1i  %10.4f %1i\n',[(1:3)' vc Xvals(xc,2) v(ii) Avals(x(ii))]')
 
% create plots for control and value
close all 

figure(1); clf
C=[0.2;0.4;0.6;0.8]*ones(1,3); colormap(C)
beliefplot(b,Avals(x),1:4);
legend({'NN','MN','NT','MT'},'location','eastoutside')

figure(2); clf
beliefplot(b,v);
colorbar
