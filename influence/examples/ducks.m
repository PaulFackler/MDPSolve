% ducks
% Based on: David R. Anderson. 1975, "Optimal Exploitation Strategies for 
%       an Animal Population in a Markovian Environment: A Theory and an
%       Example" Ecology, 56(6): 1281-1297
%
% Note that the results here differ significantly from the reported results
% Compare table with Table 1 (p. 1290) 
close all
clc
disp('Anderson (1975) duck harvest model with alternative models')

%%
mu     = 16.46;       % rainfall mean
sigma2 = 4.41;        % rainfall variance
delta  = 0.999999;    % discount factor 
% duck transition function (see ducks)
young  = @(Nn,Pn) 1./( (1./(12.48*Pn.^0.851))+(0.519./Nn) ) ;   % Eqn. 3
fall   = @(Nn,Yn) 0.92*Nn + Yn;                                 % Eqn. 5 
harv   = @(Dn,Fn) ifthenelse(Fn>0,min(Dn./Fn,1),0);
% Additive model
surv   = @(Nn,Yn,Hn) Nn.*(1-0.37*exp(2.78*Hn)) + ...
                     Yn.*(1-0.49*exp(0.90*Hn));                 % Eqns.7,8 
Nplus  = @(Nn,Pn,Dn) surv(Nn,young(Pn,Nn),harv(Dn,fall(Nn,young(Nn,Pn)))); 
% pond transition function
Pplus  = @(Pn,Rn) -2.76 + 0.391*Pn + 0.233*Rn; 
% reward function
Reward = @(N,D) ifthenelse(D<=N,D,-inf);

%%
% state and action grid definitions
nmin=3;   nmax=18;    ninc=0.5;
pmin=0.5; pmax=3.5;   pinc=0.5;
dmin=0;   dmax=12;    dinc=0.05;

nr=11;   % number of random rainfall values

%% setup model 
n=(nmin:ninc:nmax)';     % population values
p=(pmin:pinc:pmax)';     % pond number values
d=(dmin:dinc:dmax)';     % harvest values
nn=length(n); np=length(p); nd=length(d);

%%
[r,pr]=qnwnorm(nr,mu,sigma2);     % rainfall random values and probabilities
X=rectgrid(n,p,d);
R=Reward(X(:,1),X(:,3));
g=@(X,e) [Nplus(X(:,1),X(:,2),X(:,3)) Pplus(X(:,2),e)];
P=g2P(g,{n,p},X,r,pr);

%%
clear model
model.P=P;
model.reward=R;  % put penalty on harvest size above current pop.
model.discount=delta;
model.X=X;
model.svars=[1 2];

%%
results=mdpsolve(model);
Xopt=X(results.Ixopt,:);

%%
% plot optimal action
figure(1); clf
patchplot(Xopt(:,1),Xopt(:,2),Xopt(:,3));
title('Optimal Harvest Level for Additive Model')
xlabel('Population (N)')
ylabel('Pond Number (P)')
colorbar

% display optimal action in a table
disp('Optimal control for additive model  (rows=N, cols=P)')
hopta=reshape(Xopt(:,3),np,nn)';
fprintf('       ');
for j=1:np, fprintf('%3.1f   ',p(j)); end; 
fprintf('\n')
for i=1:nn
  fprintf('%4.1f   ',n(i));
  for j=1:np, fprintf('%4.2f  ',hopta(i,j)); end; 
  fprintf('\n')
end

