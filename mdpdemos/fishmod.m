clc
disp('Optimal fish harvest example with modified policy iteration')
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
 g = @(X) recruitment(X(:,1),X(:,2));

clear options
options.algorithm='f';
options.modpol=0;

clear model1
model1.name='Fish harvest';
model1.X=X;
model1.P=g2P(g,s,X);
model1.R=harvest(X(:,1),X(:,2));
model1.discount=0.99;
model1.ns=ns;
model1.Ix=Ix;
results1=mdpsolve(model1,options);


options.modpol=500;
results2=mdpsolve(model1,options);

disp('Run time and number of main iterations')
disp('ordinary function iteration')
fprintf('%12.4f  %6i\n',[results1.time  results1.iter])
disp('modified policy iteration (with maximum # of sub-iterations = 500)')
fprintf('%12.4f  %6i\n',[results2.time  results2.iter])
