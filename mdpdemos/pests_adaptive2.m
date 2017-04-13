% adaptive management of pest infestation
% 3 infestation levels - action is treat/don't treat
% 3 alternative models for non-treatment action
clc
%
close all
disp('3 Model Adaptive Management of a Pest Infestation')
p=50; % number of belief intervals

delta=0.95;   % discount factor
D=[0 10 20]'; % damage costs
Ct=10;        % treatment cost
T=inf;        % time horizon

% reward matrix (negative cost matrix)
R=-[D D+Ct];
% non-treatment transition matrix - model 1
Pn1=[0.9 0.6 0.0
     0.1 0.3 0.5
     0.0 0.1 0.5];
% non-treatment transition matrix - model 2
Pn2=[0.7 0.3   0
     0.25 0.4 0.25
     0.05 0.3 0.75];
% non-treatment transition matrix - model 3
Pn3=[0.5 0.0 0.0
     0.4 0.5 0.0
     0.1 0.5 1.0];
% treatment transition matrix
Pt=[0.9 0.8 0.6
    0.1 0.2 0.4
    0.0 0.0 0.0];

P1=[Pn1 Pt];
P2=[Pn2 Pt];
P3=[Pn3 Pt];
P={P1,P2,P3};

disp('P1'); disp(P1)
disp('P2'); disp(P2)  
disp('P3'); disp(P3) 

clear options
options.print=0;
% no model uncertainty setup
clear model
model.reward=R;
model.discount=delta;
Svals=[1;2;3];
Avals=[1;1;1;2;2;2];
Ix=[1;2;3;1;2;3];

% set up the belief state problem
[b,Pb,Rb,Svalsb,Avalsb,Ixb]=amdp(p,P,R,Svals,Avals,Ix);
clear modelA
modelA.name='Active Adaptive Management Model';
modelA.R=Rb;
modelA.P=Pb;
modelA.discount=delta;
modelA.Ix=Ixb;

% call basic solver
results=mdpsolve(modelA,options);
v=results.v; x=results.Ixopt;

% set up and solve perfect certainty models
model1=model; model1.name='Perfect Certainty that Model 1 is correct'; model1.P=P{1};
model2=model; model2.name='Perfect Certainty that Model 2 is correct'; model2.P=P{2};
model3=model; model3.name='Perfect Certainty that Model 3 is correct'; model3.P=P{3};
results1=mdpsolve(model1,options); v1=results1.v; x1=results1.Ixopt;
results2=mdpsolve(model2,options); v2=results2.v; x2=results2.Ixopt;
results3=mdpsolve(model3,options); v3=results3.v; x3=results3.Ixopt;
disp('State Values and Perfect Certainty Optimal Actions')
disp([Svals Avals(x1) Avals(x2) Avals(x3)])

ind=Svalsb(:,1)==2;
figure(1); clf
colormap([.6 .6 .6;.3 .3 .3])
set(gcf,'units','normalized','outerposition',[0.6,0.4,.4,0.5])
patchplot(Svalsb(ind,2),Svalsb(ind,3),Avalsb(x(ind),end),[1 2]);
h=legend({'don''t treat','treat'},'location','southoutside','orientation','horizontal');
pos=get(h,'position'); pos(2)=0.025; set(h,'position',pos)
pos=get(gca,'outerposition'); pos(2)=0.1; pos(4)=pos(4)*0.9; set(gca,'outerposition',pos)
title('State 2 Optimal Treatments')
xlabel('Prob(model=1)')
ylabel('Prob(model=2)')
axis square