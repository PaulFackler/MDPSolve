clc
disp('Optimal fish Harvest example with use of EV option')
c1=1.15;
c2=0.15;

recruitment = @(S,A) c1*S./(1+c2*S).*(1-A);
harvest     = @(S,A) c1*S./(1+c2*S).*A;

ns=126;
smax=1.25;
na=50;

s=linspace(0,smax,ns)';
a=linspace(0,1,na)';
X=rectgrid(s,a);
[Ix,S]=getI(X,1);
Splus=recruitment(X(:,1),X(:,2));
EV = @(V) pchip(s,V,Splus);
 g = @(X) recruitment(X(:,1),X(:,2));

clear options
options.algorithm='f';
options.modpol=500;
options.getAopt=1;

clear model1
model1.name='Fish EV demo - uses pchip interpolation to compute E[V+|X]';
model1.X=X;
model1.P=EV;
model1.EV=true;
model1.R=harvest(X(:,1),X(:,2));
model1.discount=0.99;
model1.ns=ns;
model1.Ix=Ix;
tic
results1=mdpsolve(model1,options);
toc
v1=results1.v; a1=results1.Xopt(:,2);


model2=model1;
model2.name='Fish EV demo - uses g2P to compute P matrix';
model2.P=g2P(g,s,X);
model2.EV=false;
tic
results2=mdpsolve(model2,options);
toc
v2=results2.v; a2=results2.Xopt(:,2);

fprintf('Maximum difference in value function: %1.4e\n',max(abs(v1-v2)))
fprintf('Maximum difference in strategy:       %1.4e\n',max(abs(a1-a2)))

ss=linspace(0,1.5,16)';
figure(1)
plot(ss,g([ss,zeros(size(ss,1),1)]),'-k.',ss,ss,'k')
set(gca,'Xtick',ss,'YTick',ss)
grid on
xlabel('S')
ylabel('S^*')

figure(2)
plot(s,[v1 v2])
xlabel('population size (S)')
ylabel('value (V)')

figure(3)
plot(s,[a1 a2])
xlabel('population size (S)')
ylabel('optimal harvest rate (A)')