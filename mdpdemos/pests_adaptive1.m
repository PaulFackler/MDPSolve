% adaptive management of pest infestation
% 3 infestation levels - action is treat/don't treat
% 2 alternative models for non-treatment action
clc
%
disp('2 Model Adaptive Management of a Pest Infestation')
p=100; % number of belief intervals

delta=0.95;   % discount factor
D=[0;10;20];  % damage costs
Ct=10;        % treatment cost
T=inf;        % time horizon

% reward matrix (negative cost matrix)
R=-[D D+Ct];

% non-treatment transition matrix - model 1
Pn1=[0.9 0.6 0.0
     0.1 0.3 0.5
     0.0 0.1 0.5];
  
% non-treatment transition matrix - model 2
Pn2=[0.5 0.0 0.0
     0.4 0.5 0.0
     0.1 0.5 1.0];

% treatment transition matrix
Pt=[0.9 0.8 0.6
    0.1 0.2 0.4
    0.0 0.0 0.0];

P1=[Pn1 Pt];
P2=[Pn2 Pt];

disp('P1'); disp(P1)
disp('P2'); disp(P2)

Svals=[1;2;3];
Avals=[1;1;1;2;2;2];
Ix=[1;2;3;1;2;3];

clear options
options.getAopt=false;
options.print=0;

% set up the belief state problem
[b,Pb,Rb,Svalsb,Avalsb,Ixb]=amdp(p,{P1, P2},R,Svals,Avals,Ix);

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
clear model
model.R=R(:);
model.P={P1, P2};
model.discount=delta;
model.Ix=Ix;
model1=model; model1.name='Perfect Certainty that Model 1 is correct'; model1.P=P1;
model2=model; model2.name='Perfect Certainty that Model 2 is correct'; model2.P=P2;
results1=mdpsolve(model1,options);
v1=results1.v; x1=results1.Ixopt;
results2=mdpsolve(model2,options);
v2=results2.v; x2=results2.Ixopt;
disp('States and Perfect Certainty Optimal Actions for Models 1 and 2')
disp([Svals Avals(x1) Avals(x2)])

% create control rule plot
close all
figure(1); clf
set(gcf,'units','normalized','outerposition',[0.6,0.4,.35,0.5])
colormap([.6 .6 .6;.3 .3 .3])
options=struct('addlegend',1,'vertical',0);
options.legendlabels={'0','1'};
ind=Svalsb(:,1)==2;
patchplot(Svalsb(Ixb(x),1),Svalsb(Ixb(x),2),Avalsb(x),[1 2]);
h=legend({'don''t treat','treat'},'location','southoutside','orientation','horizontal');
pos=get(h,'position'); pos(2)=0.025; set(h,'position',pos)
pos=get(gca,'outerposition'); pos(2)=0.1; pos(4)=pos(4)*0.9; set(gca,'outerposition',pos)
xlabel('Infestation Level')
ylabel('Prob(model=1)')
set(gca,'xtick',[1 2 3])