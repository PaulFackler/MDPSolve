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

delta  = 0.98;    % discount factor 

rstd=sqrt(sigma2);
% duck transition function (see ducks)
young  = @(Nn,Pn) max(0,1./( (1./(12.48*Pn.^0.851))+(0.519./Nn) )) ;   % Eqn. 3
fall   = @(Nn,Yn) max(0,0.92*Nn + Yn);                                 % Eqn. 5 
harv   = @(Dn,Fn) ifthenelse(Fn>0,min(Dn./Fn,1),0);
% Additive model
surv   = @(Nn,Yn,Hn) max(0, Nn.*(1-0.37*exp(2.78*Hn)) + ...
                           Yn.*(1-0.49*exp(0.90*Hn)));            % Eqns.7,8 
Nplus  = @(Nn,Pn,Dn) max(0,surv(Nn,young(Nn,Pn),harv(Dn,fall(Nn,young(Nn,Pn))))); 
% pond transition function
Pplus  = @(Pn,Rn) max(0,-2.76 + 0.391*Pn + 0.233*Rn); 
% reward function
Reward = @(N,D) ifthenelse(D<=N,D,-inf);

%%
% state and action grid definitions
nmin=1;   nmax=25;    ninc=0.25;
pmin=0;   pmax=3.5;   pinc=0.25;
dmin=0;   dmax=11;    dinc=0.05;

nr=31;   % number of random rainfall values

%% setup model 
n=(nmin:ninc:nmax)';     % population values
p=(pmin:pinc:pmax)';     % pond number values
d=(dmin:dinc:dmax)';     % harvest values
nn=length(n); np=length(p); nd=length(d);

%%
rprob=rvdef('ne',[mu;sqrt(sigma2)],nr);
D=[];
D=add2diagram(D,'N','s',true,{},n,[0.2198,0.7868],[]);
D=add2diagram(D,'P','s',true,{},p,[0.2198,0.5581],[]);
D=add2diagram(D,'D','a',true,{},d,[0.2107,0.4195],[]);
D=add2diagram(D,'r','c',true,{},rprob,[0.6400,0.2745],[]);
D=add2diagram(D,'Y','c',true,{'N','P'},young,[0.3994,0.6239],[5 1;5 1]);
D=add2diagram(D,'F','c',true,{'N','Y'},fall,[0.5997,0.6237],[5 1;5 1]);
D=add2diagram(D,'H','c',true,{'D','F'},harv,[0.6621,0.4175],[5 1;3 8]);
D=add2diagram(D,'N+','f',true,{'N','Y','H'},surv,[0.8298,0.7851],[5 1;5 1;6 1]);
D=add2diagram(D,'P+','f',true,{'P','r'},Pplus,[0.8350,0.5597],[5 1;5 1]);
D=add2diagram(D,'U','u',true,{'N','D'},Reward,[0.5188,0.1784],[5 1;5 1]);

%%
foptions=struct(...
'name','',...
'figpos',[0.156 0.241 0.796 0.664],...
'fontsize',0.032,...
'fontname','Rockwell',...
'backgroundcolor',[0.850 0.850 0.850],...
'nodecolor',[0.000 0.950 0.900]);
figure(1); clf
drawdiagram(D)

%%
options = struct('cleanup',0);
model=d2model(D,options);
model.d=delta;

%%
results=mdpsolve(model);
Xopt=model.X(results.Ixopt,:);

%%
% plot optimal action
figure(2); clf
patchplot(Xopt(:,2),Xopt(:,3),Xopt(:,1));
title('Optimal Harvest Level for Additive Model')
xlabel('Population (N)')
ylabel('Pond Number (P)')
colorbar

% display optimal action in a table (Compare with Table 1 in Anderson)
% results do not match
disp('Optimal control for additive model  (rows=N, cols=P)')
hopta=reshape(Xopt(:,1),np,nn)';
fprintf('       ');
for j=3:2:np, fprintf('%3.1f   ',p(j)); end; 
fprintf('\n')
for i=21:4:69
  fprintf('%4.1f   ',n(i));
  for j=3:2:np, fprintf('%4.2f  ',hopta(i,j)); end; 
  fprintf('\n')
end
%%

% longrun results
pstar = results.pstar; 
pstar(pstar<0)=0;
pstar= normalizeP(pstar);
lrp = longrunP(pstar);
m = marginals(lrp,[nn np]);

disp('longrun means of X and S')
disp(lrp'*Xopt)
disp([m{1}'*n, m{2}'*p])
AA=Xopt(:,1); 
s0=lrp'*Xopt(:,2:3);

optharv = @(S) rectbas(S,{n,p},[],1)'*Xopt(:,1);
noharvest = @(S) zeros(size(S,1),1);

T=50;
s0=[8 2];
reps=10000;
s0=repmat(s0,reps,1);
[Y1,z]=dsim(D,s0,T+1,optharv);
Y2=dsim(D,s0,T+1,noharvest,[],z);
figure(3); clf;
set(gcf,'units','normalized','position',[0.219 0.231 0.620 0.533])
plot(0:T,[mean(Y2{1},1);mean(Y1{1},1)],0:T,mean(Y1{3},1),'LineWidth',2)
xlabel('t')
legend('Mean population w/ no harvest','Mean population w/ optimal harvest','Mean Optimal harvest','location','eastoutside')


%%
[probs,pparents] = gs2ps({Pplus,Nplus},{[1 -1],[2 1 3]},{p,n},{p,n,d},{rprob},options);
EV=EVcreate(probs,pparents,{p,n,d});
X=rectgrid({p,n,d});
R = Reward(X(:,2),X(:,3));
Ix = getI(X,[1 2]);
model = struct('P',EV,'R',R,'d',delta,'Ix',Ix,'X',X,'EV',1);

moptions=struct('algorithm','i');
results=mdpsolve(model,moptions);
Xopt=model.X(results.Ixopt,:);


figure(4); clf
patchplot(Xopt(:,2),Xopt(:,1),Xopt(:,3));
title('Optimal Harvest Level for Additive Model')
xlabel('Population (N)')
ylabel('Pond Number (P)')
colorbar
