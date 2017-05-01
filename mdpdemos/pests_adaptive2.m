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
Pn2=[0.7  0.3   0
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

clear model
model.reward=R;
model.discount=delta;

Svals=[1;2;3];
Avals=[1;2];
X=rectgrid(Avals,Svals);
Ix=getI(X,2);

clear options
options.print=0;

%% active adaptive management
% set up the belief state problem
[b,Pb,Rb,Svalsb,Xb,Ixb]=amdp(p,{P1,P2,P3},R(:),Svals,X,Ix);

clear modelA
modelA.name='Active Adaptive Management Model';
modelA.R=Rb;
modelA.P=Pb;
modelA.discount=delta;
modelA.Ix=Ixb;

% call basic solver
results=mdpsolve(modelA,options);
v=results.v; Xopt=Xb(results.Ixopt,:);

% set up and solve perfect certainty models
model1=model; model1.name='Perfect Certainty that Model 1 is correct'; model1.P=P{1};
model2=model; model2.name='Perfect Certainty that Model 2 is correct'; model2.P=P{2};
model3=model; model3.name='Perfect Certainty that Model 3 is correct'; model3.P=P{3};
results1=mdpsolve(model1,options); v1=results1.v; x1=results1.Ixopt;
results2=mdpsolve(model2,options); v2=results2.v; x2=results2.Ixopt;
results3=mdpsolve(model3,options); v3=results3.v; x3=results3.Ixopt;
disp('State Values and Perfect Certainty Optimal Actions')
disp([Svals X(x1,1) X(x2,1) X(x3,1)])

ind=Xopt(:,2)==2;
figure(1); clf; 
set(gcf,'units','normalized','position',[0.5    0.5    0.35    0.4])
beliefplot(Xopt(ind,3:5),Xopt(ind,1),[1 2]); 
pos=get(gca,'position');
pos(2)=pos(2)-.04; 
set(gca,'position',pos)
h=title('Optimal Action for S=2');
pos=get(h,'position');
pos(2)=pos(2)+.025; 
set(h,'position',pos,'fontsize',12)
legend('Do nothing','Treat')

figure(2); clf; 
set(gcf,'units','normalized','position',[0.25    0.55    0.65    0.35])
rv=[min(v) max(v)];
for i=1:3
  ind=Xopt(:,2)==2;
  subplot(1,3,i)
  beliefplot(Xopt(ind,3:5),v(ind),rv); 
  title(['Value Function, S=' num2str(i)])
end
h=colorbar;
pos=get(h,'position');
pos(1)=1-2.5*pos(3); pos(2)=0.25; pos(4)=0.5;
set(h,'position',pos)