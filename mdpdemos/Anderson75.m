% ducks
% Based on: David R. Anderson. 1975, "Optimal Exploitation Strategies for 
%       an Animal Population in a Markovian Environment: A Theory and an
%       Example" Ecology, 56(6): 1281-1297
%
% Note that the results here differ significantly from the reported results
% Compare table with Table 1 (p. 1290) (or Table 2 (p.1292) if the
% compensatory option is used).
close all
clear variables
%clc
disp('Anderson (1975) duck harvest model with alternative models')
mdpoptions=struct('print',0,'algorithm','p','vanish',0.999999,'relval',0);

delta=1;    % discount factor 
% state and action grid definitions
nmin=3;   nmax=18;    ninc=0.1;
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
hshock=[1;1];
hprob=[.5;.5];
% combine two noise terms
[e,w]=mergestates({hshock,hprob},{r,pr});

% preliminary plots
figure(1)
p=linspace(pmin,pmax,1001)';
pplus1=Pplus(p,mu-sqrt(sigma2));
pplus2=Pplus(p,mu);
pplus3=Pplus(p,mu+sqrt(sigma2));
plot(p,[pplus1 pplus2 pplus3],p,p,':k')
xlabel('Number of Ponds (P)')
ylabel('P^+')
title('Pond Number Transition Function')
legend({'R=\mu-\sigma','R=\mu','R=\mu+\sigma'},'location','southeast')

p=(pmin:pinc:pmax)';
p=linspace(pmin,pmax,51)';
P=g2P(Pplus,p,p,r,pr);
P(P<0)=0; P=mxv(P,1./sum(P));
lrP=longrunP(P);

a=-2.76; b=0.391; c=0.233; 
Ep=(a+c*mu)/(1-b);
Sp=c*sqrt(sigma2/(1-b^2));
lrP2=pdfn(p,Ep,Sp); lrP2=lrP2/sum(lrP2);

figure(2)
plot(p,lrP,'-k.',p,lrP2)
title('Correct and Model Based Long Run Pond Number Distribution')
xlabel('Pond Number (P)')

figure(3)
n=linspace(eps,20,1001)';
nt=@(n,h) Anderson75tran(n,Ep,h,1,0);
plot(n,[nt(n,0) nt(n,2) nt(n,4) nt(n,6) nt(n,8)],n,n,':k')
xlabel('Population Size (N)')
ylabel('N^+')
title('Compensatory Population Size Transition Function (with average pond numbers)')
legend({'h=0','h=2','h=4','h=6','h=8'},'location','southeast')

figure(4)
nt=@(n,h) Anderson75tran(n,Ep,h,1,1);
plot(n,[nt(n,0) nt(n,2) nt(n,4) nt(n,6) nt(n,8)],n,n,':k')
xlabel('Population Size (N)')
ylabel('N^+')
title('Additive Population Size Transition Function (with average pond numbers)')
legend({'h=0','h=2','h=4','h=6','h=8'},'location','southeast')

% Solve Model 
n=(nmin:ninc:nmax)';     % population values
p=(pmin:pinc:pmax)';     % pond number values
d=(dmin:dinc:dmax)';     % harvest values
nn=length(n); np=length(p); nd=length(d);

X=rectgrid(n,p,d);
N=X(:,1); P=X(:,2); D=X(:,3);
[Ix,S]=getI(X,[1 2]);
clear model
model.reward=reshape(D-(N<=D & D>0)*1e10,nn*np,nd);  % put penalty on harvest size above current pop.
model.discount=delta;
model.Ix=Ix;
model.colstoch=1;
model.n=size(S,1);
g2Poptions=struct('cleanup',0);

%% compensatory model
transfunc=@(x,e) [Anderson75tran(x(:,1),x(:,2),x(:,3),e(:,1),0) Pplus(x(:,2),e(:,2))];
model.name='Anderson (1975) - compensatory model';
model.P=g2P(transfunc,{n,p},X,e,w,g2Poptions);
tic
results=mdpsolve(model,mdpoptions);
toc
vc=results.v; xc=results.Ixopt; pstarc=results.pstar;

% additive model
transfunc=@(x,e) [Anderson75tran(x(:,1),x(:,2),x(:,3),e(:,1),1) Pplus(x(:,2),e(:,2))];
model.name='Anderson (1975) - additive model';
model.P=[];
model.P=g2P(transfunc,{n,p},X,e,w,g2Poptions);
tic
results=mdpsolve(model,mdpoptions);
toc
va=results.v; xa=results.Ixopt; pstara=results.pstar;
%%

options=struct(...
      'clim', [0 9], ...
      'grayscale',    0, ...
      'squareplot',   0, ...
      'addlegend',    1, ...
      'vertical',     1, ...     % 1 if legend is vertical
      'colorbartype', 1);       

clim=[0 max(max(D(xc)),max(D(xa)))];    
figure(5); clf
mdpplot(S,D(xc),[1 2],{'N','P'},options);
title('Optimal Control for Compenstory Model')

figure(6); clf
mdpplot(S,D(xa),[1 2],{'N','P'},options);
title('Optimal Control for Additive Model')

disp('Optimal control for compensatory model (rows=N, cols=P)')
hoptc=reshape(D(xc),np,nn)'; 
  fprintf('       ');for j=1:np, fprintf('%3.1f   ',p(j)); end; fprintf('\n')
for i=1:nn
  fprintf('%4.1f   ',n(i));for j=1:np, fprintf('%4.2f  ',hoptc(i,j)); end; fprintf('\n')
end

disp('Optimal control for additive model  (rows=N, cols=P)')
hopta=reshape(D(xa),np,nn)';
  fprintf('       ');for j=1:np, fprintf('%3.1f   ',p(j)); end; fprintf('\n')
for i=1:nn
  fprintf('%4.1f   ',n(i));for j=1:np, fprintf('%4.2f  ',hopta(i,j)); end; fprintf('\n')
end
%%
figure(7)
plot(n,hoptc(:,[find(p==0.5) find(p==1.5) find(p==2.5) find(p==3.5)]))
xlabel('N')
ylabel('D')
legend({'0.5','1.5','2.5','3.5'},'location','northwest')
ylim([0 9])
title('Optimal Harvest at Alternative Pond Levels - Compensatory Model')

figure(8)
plot(n,hopta(:,[find(p==0.5) find(p==1.5) find(p==2.5) find(p==3.5)]))
xlabel('N')
ylabel('D')
legend({'0.5','1.5','2.5','3.5'},'location','northwest')
ylim([0 9])
title('Optimal Harvest at Alternative Pond Levels - Additive Model')



% Long-run results
plrc=longrunP(pstarc);  % long run state probabilitiy distribution (compensatory)
plra=longrunP(pstara);  % long run state probabilitiy distribution (additive) 
disp('Long run average harvest size (using E[v] and (1-delta)*E[v]')
disp('Compensatory')
disp([plrc'*D(xc), plrc'*vc*(1-delta)])
disp('Additive')
disp([plra'*D(xa), plra'*va*(1-delta)])

Mc=marginals(plrc,[nn np]);
p1c=Mc{1}; p2c=Mc{2};
Ma=marginals(plra,[nn np]);
p1a=Ma{1}; p2a=Ma{2};
 
figure(7)
lrP2=pdfn(p,Ep,Sp); lrP2=lrP2/sum(lrP2);
plot(p,[p2c p2a],'.-',p,lrP2)
title('Model Based Long Run Pond Number Distribution')
xlabel('Pond Number (P)')
 
figure(8)
plot(n,[p1c p1a],'.-')
title('Long run distribution of Population Size')
xlabel('N')
legend({'Compensatory','Additive'},'location','northeast')

disp('Longrun mean population level for compensatory and additive models')
disp([p1c'*n  p1a'*n])

disp('Estimate of Long-run Population Std. Dev. for compensatory and additive models')
disp([sqrt(p1c'*(n-p1c'*n).^2)  sqrt(p1a'*(n-p1a'*n).^2)])
