% ducks
% Based on: David R. Anderson. 1975, "Optimal Exploitation Strategies for 
%       an Animal Population in a Markovian Environment: A Theory and an
%       Example" Ecology, 56(6): 1281-1297
%
% Note that the results here differ significantly from the reported results
% Compare table with Table 1 (p. 1290) (or Table 2 (p.1292) if the
% compensatory option is used).
clc
close all
%clear variables

disp('Anderson (1975) duck harvest model with alternative models')
mdpoptions=struct('print',2,'algorithm','p','vanish',0.999999,'relval',0);

delta=.98;    % discount factor 
% state and action grid definitions
nmin=2;   nmax=18;    ninc=0.1;
pmin=0.5; pmax=3.5;   pinc=0.1;
dmin=0;   dmax=12;    dinc=0.05;

nr=11;   % number of random rainfall values for Gaussian quadrature

Pplus = @(Pn,Rn) min(pmax,max(pmin,-2.76 + 0.391*Pn + 0.233*Rn));    % pond transition function
Pplus = @(Pn,Rn) -2.76 + 0.391*Pn + 0.233*Rn;    % pond transition function

mu=16.46;        % rainfall mean
sigma2=4.41;     % rainfall variance
errorcase=2;     % (0) use expected value 
                 % (1) use Gaussian quadrature 
                 % (2) use Lubow's discretization
switch errorcase
  case 0
    r=mu; pr=1;
  case 1
   [r,pr]=qnwnorm(nr,mu,sigma2);     % rainfall random values and probabilities
   %Ep=(a+c*mu)/(1-b);                % exact long run mean pond number
  case 2
    % alternative rain distribution from Lubow's ASDP User's Guide
    r=[11.58;13.01;13.77;14.29;14.69;15.05;15.36;15.65;15.93;16.20;16.46; ...
       16.72;16.99;17.27;17.56;17.87;18.23;18.63;19.15;19.91;21.34];
    pr=ones(21,1)*0.05; pr(1)=0.025; pr(end)=0.025;
end
nr=length(r);
% harvest shock
hshock=[0.9;1;1.1];
hprob=[0.2;0.6;0.2];
% combine two noise terms
[e,w]=mergestates({hshock,hprob},{r,pr});


% Solve Model 
n=(nmin:ninc:nmax)';     % population values
p=(pmin:pinc:pmax)';     % pond number values
d=(dmin:dinc:dmax)';     % harvest values
nn=length(n); np=length(p); nd=length(d);

X=rectgrid(n,p,d);
feasible=X(:,3)<=X(:,1);
X=X(feasible,:);
[Ix,S]=getI(X,[1 2]);


%% compensatory model
transfunc=@(x,e) [Anderson75tran(x(:,1),x(:,2),x(:,3),e(:,1),0) Pplus(x(:,2),e(:,2))];

clear model
model.name='Anderson (1975) - compensatory model';
goptions=struct('cleanup',1);
tic
model.P=g2P(transfunc,{n,p},X,e,w,goptions);
toc
%%
model.reward=X(:,3);  
model.discount=delta;
model.Ix=Ix;
model.colstoch=1;
model.n=size(S,1);

mdpoptions=struct('print',2);
tic
results=mdpsolve(model,mdpoptions);
toc
v=results.v; Xopt=X(results.Ixopt,:); 

figure(1); clf
patchplot(Xopt(:,1),Xopt(:,2),Xopt(:,3));
title('Optimal Control for Compensatory Model')
