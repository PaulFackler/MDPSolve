clearclose
disp('Anderson (1975) duck harvest model using g2EV')

% model parameters
delta  = 0.98;               % discount factor 
mu     = 16.46;              % rainfall mean
sigma2 = 4.41;               % rainfall variance
% pond transition function
Pplus = @(Pn,Rn) -2.76 + 0.391*Pn + 0.233*Rn;  
% population transition function
Nplus = @(N,P,D,Hn) Anderson75tran(N,P,D,Hn,0);

% harvest shock
hshock=[0.9;1;1.1];
hprob=[0.2;0.6;0.2];

% solution parameters
nmin=2;   nmax=18;    ninc=0.1;
pmin=0.5; pmax=3.5;   pinc=0.1;
dmin=0;   dmax=12;    dinc=0.05;
nr=11;          % number of random rainfall noise values
cleanup=2;      % project to boundary

% get variable values and indices
rrv=rvdef('n',[mu;sqrt(sigma2)],nr);
hrv=rvdef('d',hprob,hshock);
e={rrv,hrv};

n=(nmin:ninc:nmax)';     % population values
p=(pmin:pinc:pmax)';     % pond number values
d=(dmin:dinc:dmax)';     % harvest values
nn=length(n); np=length(p); nd=length(d);

X=rectgrid(n,p,d);
feasible=X(:,3)<=X(:,1);
X=X(feasible,:);
[Ix,S]=getI(X,[1 2]);

% use g2EV to get EV function
goptions=struct('cleanup',cleanup,'order',[2 1]);
%goptions=struct('cleanup',cleanup);
xelist={[1 2 3 -2],[2 -1]};
tic
[EV,pp,Is,ws]=g2EV({Nplus Pplus},{n,p},X,e,xelist,goptions);
fprintf('time to obtain EV function:       %8.4f\n',toc)
% demonstrates how P can be obtained using p
tic
ii=getI(X,2); PP=kroncol(pp{1},pp{2}(:,ii));  
fprintf('time to obtain P from EV results: %8.4f\n',toc)

% use g2P to get transition matrix
transfunc=@(x,e) [Anderson75tran(x(:,1),x(:,2),x(:,3),e(:,2),0) Pplus(x(:,2),e(:,1))];
goptions=struct('cleanup',cleanup);
tic
P=g2P(transfunc,{n,p},X,e,[],goptions);
fprintf('time to obtain P matrix:          %8.4f\n',toc)

disp('check size of difference between P matrices obtained in two different ways')
disp(full(max(max(abs(P-PP)))))
%%
clear model
model.P=EV;
model.EV=1;
model.reward=X(:,3);  
model.discount=delta;
model.Ix=Ix;
model.colstoch=1;
model.n=size(S,1);
mdpoptions=struct('print',0,'algorithm','i','debug',1);

% use EV function to solve model
tic
results=mdpsolve(model,mdpoptions);
fprintf('time used by mdpsolve:            %8.4f\n',toc)
Xopt1=X(results.Ixopt,:); v1=results.v;

% use P matrix to solve model
model.P=P;
model.EV=0;
tic
results=mdpsolve(model,mdpoptions);
fprintf('time used by mdpsolve:            %8.4f\n',toc)
Xopt2=X(results.Ixopt,:);  v2=results.v;

%%
figure(1); clf
colormap(linspace(.95,.05,64)'*[1 1 1]);
set(gcf,'units','normalized','position',[0.203 0.241 0.681 0.464])
subplot(1,2,1)
patchplot(Xopt2(:,1),Xopt2(:,2),Xopt2(:,3)-Xopt1(:,3));
title('difference in action')
xlabel('population (N)')
ylabel('ponds (P)')
colorbar
subplot(1,2,2)
patchplot(Xopt2(:,1),Xopt2(:,2),v1-v2);
title('difference in value')
xlabel('population (N)')
ylabel('ponds (P)')
colorbar