% adaptive management of pest infestation
% 3 infestation levels - action is treat/don't treat
% 2 alternative models for non-treatment action
clc
%
disp('2 Model Adaptive Management of a Pest Infestation')
p=100; % number of belief intervals

delta=0.95;        % discount factor
Damage=[0;10;20];  % damage costs
Ct=10;             % treatment cost
T=inf;             % time horizon

% reward matrix (negative cost matrix)
R=-[Damage Damage+Ct];

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
Avals=[1;2];
X=rectgrid(Avals,Svals);
Ix=getI(X,2);

clear options
options.print=0;

%% active adaptive management
% set up the belief state problem
[b,Pb,Rb,Svalsb,Xb,Ixb]=amdp(p,{P1, P2},R(:),Svals,X,Ix);

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
disp([Svals X(x1,1) X(x2,1)])

% create control rule plot
close all
figure(1); clf
set(gcf,'units','normalized','outerposition',[0.3,0.4,.7,0.5])
colormap([.6 .6 .6;.3 .3 .3])
subplot(1,2,1)
patchplot(Xopt(:,2),Xopt(:,3),Xopt(:,1),[1 2]);
h=legend({'don''t treat','treat'},'location','southoutside','orientation','horizontal');
pos=get(h,'position'); pos(2)=0.025; set(h,'position',pos)
pos=get(gca,'outerposition'); pos(2)=0.1; pos(4)=pos(4)*0.9; set(gca,'outerposition',pos)
xlabel('Infestation Level')
ylabel('Prob(model=1)')
set(gca,'xtick',[1 2 3])
title('Active Adaptive Strategy')

%% passive adaptive strategy (loops over b)
modelp=model;
Ap=zeros(size(b,1),3);
vp=zeros(size(b,1),3);
for i=1:size(b,1)
  bb=b(i,1);
  modelp.P=bb*P1+(1-bb)*P2;
  resultsp=mdpsolve(modelp,options);
  vp(i,:)=resultsp.v'; Ap(i,:)=X(resultsp.Ixopt,1)';
end
vp=vp(:);
Ap=Ap(:);

subplot(1,2,2)
patchplot(Xopt(:,2),Xopt(:,3),Ap,[1 2]);
h=legend({'don''t treat','treat'},'location','southoutside','orientation','horizontal');
pos=get(h,'position'); pos(2)=0.025; set(h,'position',pos)
pos=get(gca,'outerposition'); pos(2)=0.1; pos(4)=pos(4)*0.9; set(gca,'outerposition',pos)
xlabel('Infestation Level')
ylabel('Prob(model=1)')
set(gca,'xtick',[1 2 3])
title('Passive Adaptive Strategy')

%% expected value of perfect information
EVPI=zeros(size(b,1),3);
for ss=1:3
  EVPI(:,ss)=b(:,1)*v1(ss)+b(:,2)*v2(ss)-v(Svalsb(:,1)==ss);
end
figure(2); ss=3; plot(b(:,1),EVPI)
yy=get(gca,'ylim'); yy(1)=0; set(gca,'ylim',yy)
legend('s=1','s=2','s=3')
xlabel('b_1')
ylabel('EVPI')

%% uses the extended POMDP approach
X2=rectgrid([1;2],Avals,Svals); 
Z=rectgrid([1;2],Svals);
[Pb2,Rb2,Sb2,Xb2,Ixb2]=xpomdp(p,blkdiag(P1,P2),[R(:) R(:)],X2,[3],1,Z,2,1);
modelx=struct('P',Pb2,'R',Rb2,'X',Xb2,'Ix',Ixb2,'d',delta);
resultsx=mdpsolve(modelx,options);
vx=resultsx.v; Xoptx=Xb2(resultsx.Ixopt,:);

figure(3); clf
set(gcf,'units','normalized','outerposition',[0.3,0.4,.7,0.5])
colormap([.6 .6 .6;.3 .3 .3])
subplot(1,2,1)
patchplot(Xoptx(:,2),Xoptx(:,3),Xoptx(:,1),[1 2]);
h=legend({'don''t treat','treat'},'location','southoutside','orientation','horizontal');
pos=get(h,'position'); pos(2)=0.025; set(h,'position',pos)
pos=get(gca,'outerposition'); pos(2)=0.1; pos(4)=pos(4)*0.9; set(gca,'outerposition',pos)
xlabel('Infestation Level')
ylabel('Prob(model=1)')
set(gca,'xtick',[1 2 3])
title('Active Adaptive Strategy')

%% diagramatic approach
Mvals=[1;2];
pesttran=rvdef('d',[P1 P2],Svals);
U=@(S,A) -Damage(S)-Ct*(A==2);

D=[];
D=add2diagram(D,'peststate','s',true,{},Svals);
D=add2diagram(D,'treat','a',true,{},Avals);
D=add2diagram(D,'model','p',false,{},Mvals);
D=add2diagram(D,'peststate+','f',true,{'model','treat','peststate'},pesttran);
D=add2diagram(D,'utility','u',true,{'peststate','treat'},U);
D.locs=[ ...
0.205 0.205 0.309 0.612 0.590;
0.718 0.551 0.383 0.718 0.478]';

doptions=struct('inc',p,'d',delta,'reps',0000);
modelD=d2model(D,doptions);

[resultsD]=mdpsolve(modelD,options);
XoptD=modelD.X(resultsD.Ixopt,:);

% create control rule plot
figure(4); clf
set(gcf,'units','normalized','outerposition',[0.3,0.4,.7,0.5])
colormap([.6 .6 .6;.3 .3 .3])
subplot(1,2,1)
patchplot(XoptD(:,2),XoptD(:,3),XoptD(:,1),[1 2]);
h=legend({'don''t treat','treat'},'location','southoutside','orientation','horizontal');
pos=get(h,'position'); pos(2)=0.025; set(h,'position',pos)
pos=get(gca,'outerposition'); pos(2)=0.1; pos(4)=pos(4)*0.9; set(gca,'outerposition',pos)
xlabel('Infestation Level')
ylabel('Prob(model=1)')
set(gca,'xtick',[1 2 3])
title('Active Adaptive Strategy')

figure(5)
drawdiagram(D)
