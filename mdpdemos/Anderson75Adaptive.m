% ducks
% Based on: David R. Anderson. 1975, "Optimal Exploitation Strategies for 
%       an Animal Population in a Markovian Environment: A Theory and an
%       Example" Ecology, 56(6): 1281-1297
%
% Note that the results here differ significantly from the reported results
% Compare table with Table 1 (p. 1290) (or Table 2 (p,1292) if the
% compensatory option is used).

%close all
clc
disp('Anderson (1975) duck harvest model - active adaptive')

mdpoptions=struct('print',2,'algorithm','p','relval',1,'vanish',0.999999,'modpol',500);
g2Poptions=struct('cleanup',0);

Pplus = @(Pn,Rn) -2.76 + 0.391*Pn + 0.233*Rn;    % pond transition function
delta=0.95;

% state and action grid definitions
nmin=3;   nmax=18;    ninc=.5;
pmin=0.5; pmax=3.5;   pinc=0.5;
dmin=0;   dmax=12;    dinc=0.1;
nr=21;   % number of random rainfall values if Gaussian quadrature used

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

% Solve Model 
n=(nmin:ninc:nmax)';     % population values
p=(pmin:pinc:pmax)';     % pond number values
d=(dmin:dinc:dmax)';     % harvest values
nn=length(n); np=length(p); nd=length(d);

X=rectgrid(n,p,d);
ind=X(:,1)>=X(:,3);   % eliminate state/action combinations with harvest size above current pop.
X=X(ind,:);
N=X(:,1); P=X(:,2); D=X(:,3);
[Ix,S]=getI(X,[1 2]);
R=min(N,D)*0.6+min(N,D*1.1)*0.2+min(N,D*0.9)*0.2;   %  ??? following Lubow in ASDP User's Guide

% compensatory model
transfunc=@(x,r) [Anderson75tran(x(:,1),x(:,2),x(:,3),1,0) Pplus(x(:,2),r)];
Pc=g2P(transfunc,{n,p},X,r,pr,g2Poptions);
% additive model
transfunc=@(x,r) [Anderson75tran(x(:,1),x(:,2),x(:,3),1,1) Pplus(x(:,2),r)];
Pa=g2P(transfunc,{n,p},X,r,pr,g2Poptions);

% set up adaptive problem
nb=5;
[b,Pb,Rb,Sb,Xb,Ixb]=amdp(nb,{Pc,Pa},R,S,X,Ix);

clear model
model.P =Pb;
model.reward=Rb;
model.discount=delta;
model.Ix=Ixb;
model.n=size(Pb,1);
model.transposed=1;
model.T=inf;

results=mdpsolve(model,mdpoptions);
if isempty(results.errors)
  v=results.v; x=results.Ixopt; pstar=results.pstar;
else
  return
end

options=struct(...
      'figuretitle',  'Optimal Control for Adaptive Model', ...
      'legendtitle',   'Harvest Level', ...
      'grayscale',    0, ...
      'squareplot',   0, ...
      'addlegend',    1, ...
      'vertical',     1, ... 
      'colorbartype', 1);       

figure(10); clf
h=mdpplot(Sb,Xb(x,3),[1 2 0 4],{'N','P','XXX','B'},options);
p1=get(h(1),'Position');
p2=get(h(2),'Position');

p1(2)=p1(2)-(p2(4)-p1(4));
p1(4)=p2(4);
set(h(1),'Position',p1)

bvals=unique(Sb(:,4));
hopt=Xb(x,3); 
disp('Mallard Management Adaptive Control')
disp('Optimal control (rows=N, cols=P)')

for i=1:nb+1
  ii=Sb(:,4)==bvals(i);
  hopti=reshape(hopt(ii),np,nn)'; 
  disp(['Belief in additive model: ' num2str(bvals(i))])
  disp([0 p';n hopti])
end


%ii=Xb(x,1)==6  & Xb(x,2)==2; disp([(0:0.2:1)' flipud(hopt(ii))])
%ii=Xb(x,1)==18 & Xb(x,2)==2; disp([(0:0.2:1)' flipud(hopt(ii))])
