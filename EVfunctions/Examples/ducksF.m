% Based on:
% U.S. Fish and Wildlife Service
% Adaptive Harvest Management, 2012 Hunting Season
% http://www.fws.gov/migratorybirds/NewReportsPublications/AHM/Year2012/AHMReport2012.pdf

% N: population level
% P: pond numbers
% HI: harvest intensity (1-4)
% HR: harvest rate
% HN: harvest noise
% H:  harvest
% AD: adult ducks
% JD: juvenile ducks
% K: kill rate (HR/(1-c))

%close all
%clear variables

phi   = 0.897;  % ratio F to M summer survival (p. 33)
m     = 0.5246; % proportion of males in breeding population (p. 33)
Ngoal = 8.5;    % population goal (p. 17)
CSC   = 4.75;   % closed season constraint (p. 15) [not used]
c     = 0.2;    % crippling-loss rate (p. 34)
evar  = 1.2567; % pond shock variance (std. dev. = 1.1210) (p. 34)
nvar  = 0.0280; % population variance (p. 35)
 
% survival rates (p. 34)
s0ma = 0.7896; % survival in absence of harvesting (male/additive)
s0fa = 0.6886; % survival in absence of harvesting +(female/additive)
s0mc = 0.6467; % survival in absence of harvesting (male/compensatory)
s0fc = 0.5965; % survival in absence of harvesting (female/compensatory)

% correction factors (p. 33)
gams = 0.9407; % survival correction factor 
gamr = 0.8647; % reproduction correction factor
 
% mean and standard deviation of harvest rate under alternative actions
% Table 2, p.16
mk=[0.0088;0.0552;0.1011;0.1169];
sd=[0.0019;0.0129;0.0214;0.0190];
% differential vulnerability rate relative to adult males (p. 16)
afvr=0.7191;  % adult females
jmvr=1.5407;  % juvenile males
jfvr=1.1175;  % juvenile females

strong   = 0;  % 0/1 weak/strong density dependence
additive = 1;  % 0/1 compensatory/additive survival

% pond transition function  (p. 34)
Ptran = @(P,rain) max(0,2.2127 + 0.3420*P + rain);

% reproduction functions (p. 34)
if strong  
  Repr = @(N,P) 1.1390 + 0.1376*P - 0.1131*N; % strong density dependence
else
  Repr = @(N,P) 0.7166 + 0.1083*P - 0.0373*N; %   weak density dependence
end

HR=@(HI,HN) HN.*sd(HI) + mk(HI);   % harvest rate
% survival functions
if additive
  Surv = @(K,sa) sa*max(0,(1-K));                            % additive survival
else
  Surv = @(K,sc) sc*(K<(1-sc)) + max(0,(1-K)).*(K>=(1-sc));  % compensatory survival
end

% duck transition function                
Nnext = @(N,J,h,PN) gams*N.*(m*Surv(h,s0ma) + (1-m)*(Surv(afvr*h,s0fa) + ...
                gamr*J.*(Surv(jfvr*h,s0fa) + Surv(jmvr*h,s0ma)*phi))).*exp(PN);
% function for the harvest level             
Harvest = @(N,J,HI) gams*(m/phi+(1-m)*(afvr + (gamr*(jfvr+jmvr))*J)).*mk(HI).*N;
            
Utility=@(H,N) H.*min(1,N/Ngoal); 
%Utility=@(H,N) H; 

% discretization for state variables
Nmin = 1;  Nmax = 35; Ninc = .1;
Pmin = 0;  Pmax = 9;  Pinc = .1;

% quadrature is used with the following number of nodes
nk=21; % number of harvest noise values
nn=21; % number of population noise values
np=21; % number of rain values 

N=(Nmin:Ninc:Nmax)';     % population in millions
P=(Pmin:Pinc:Pmax)';     % pond numbers in millions
HI=(1:4)';               % harvest intensity

% use RV structures for the noise terms
cptH = rvdef('n',[0;1]         ,nk); % harvest survival noise
cptN = rvdef('n',[0;sqrt(nvar)],nn); % population noise
cptR = rvdef('n',[0;sqrt(evar)],np); % rain noise

Ntran = @(HI,P,N,HN,PN) Nnext(N,Repr(N,P),HR(HI,HN),PN);
Rfunc=@(HI,P,N) Utility( Harvest(N, Repr(N,P), HI), N );
s={P,N};
x={HI,P,N};
e={cptH,cptN,cptR};
g={Ptran,Ntran};
gparents={[2 -3],[1 2 3 -1 -2]};

X=rectgrid(x);
R=Rfunc(X(:,1),X(:,2),X(:,3));
svars = [2 3];
[Ix,S] = getI(X, svars);
ns=size(S,1);
nx=size(X,1);

% Convert functions to CPTs with gs2ps and get EV function with EVcreate
options=struct('cleanup',2);
tic
[p,parents]=gs2ps(g,gparents,s,x,e,options);
fprintf('time taken to convert functions to CPTs: %1.4f\n',toc)
tic
EV=EVcreate(p,parents,x,e);
fprintf('time taken to create EV function: %1.4f\n',toc)

% set up DP model and solve
clear model
model.P  = EV;
model.EV = true;
model.R  = R;
model.Ix = Ix;
model.d  = 0.98;
options.algorithm='i'; options.print=2;
results=mdpsolve(model,options);

%%%%  alternative approaches
if 0   % construct P matrix with fm2P
  tic
  P=fm2P(p,parents,x,e);  
  fprintf('time taken to create transition matrix: %1.4f\n',toc)
  modelP=model;
  modelP.P = P;
  modelP.EV = false;
  resultsP=mdpsolve(modelP,options);
  if max(abs(results.Ixopt-resultsP.Ixopt))~=0
    disp('factor EV and diagram give different results')
  end
end

if 0   % convert to diagram with f2d and use d2model
  tic
  options.ptype=0; % create EV function
  D=f2d(g,gparents,s,x,e);
  modelDD=d2model(D,options);
  fprintf('time taken to create and process diagram: %1.4f\n',toc)
  modelDD.d  = 0.98;
  modelDD.R  = R;
  modelDD.Ix = Ix;
  resultsDD=mdpsolve(modelDD,options);
  if max(abs(results.Ixopt-resultsDD.Ixopt))~=0
    disp('factor EV and diagram give different results')
  end
end



%%
if isfield(results,'Ixopt') && ~isempty(results.Ixopt)
  Xopt=X(results.Ixopt,:);
  figure(2); clf
  set(gcf,'units','inches','position',[7  5  8 4.5])
  C=[0.2;0.4;0.6;0.8]*[1 1 1]; colormap(C)
  patchplot(Xopt(:,2),Xopt(:,3),Xopt(:,1),1:4);
  xlabel('ponds')
  ylabel('adult ducks')
  title('Optimal harvest regime strategy')
  legend('restricted','conservative','moderate','liberal','location','eastoutside')
  ylim([N(1) 12])
else
  mdpreport(results)
end

modelF=model;
return
%% plot long run marginal and joint distributions
lrp=longrunP(results.pstar,struct('fast',1)); 
M=marginals(lrp,D.sizes(1:2));  
figure(3); clf
set(gcf,'units','normalized','position',[0.15 0.2 0.8 0.45])
subplot(1,3,1);  
vname='ponds';
plot(dvalues(D,vname,'m'),M{1},'k'); 
xlabel(vname)
xlim([0 Pmax])
subplot(1,3,2);  
vname='adult ducks';
plot(dvalues(D,vname,'m'),M{2},'k'); 
xlabel(vname)
xlim([0 Nmax])
title('Marginal and joint long run probability distributions','Fontsize',16)

S=dvalues(D,1:2,'m'); 
subplot(1,3,3); 
patchplot(S(:,1),S(:,2),lrp);
xlabel('ponds')
ylabel('adult ducks')
h=colorbar;
pos=get(h,'position'); pos(1)=pos(1)+(1-pos(1))/2; set(h,'position',pos)

%% Q-functions
Q=zeros(size(Xopt,1),4);
for i=1:4
   ind = find(X(:,1)==i);
   Q(:,i) = R(ind)+model.d*EV(results.v,ind);
end

%%
  figure(4); clf
  set(gcf,'units','inches','position',[7  5  8 4.5])
  %C=[0.2;0.4;0.6;0.8]*[1 1 1]; colormap(C)
  for i=2:4
    subplot(1,3,i-1);
    patchplot(Xopt(:,2),Xopt(:,3),Q(:,i)./Q(:,1)-1,[-0.02 0.06]);
    xlabel('ponds')
    ylabel('adult ducks')
    title(['Q(S,' num2str(i) ')/Q(S,1)-1'])
    colorbar
  end
  
  %%
  figure(3); clf
  set(gcf,'units','inches','position',[7  5  8 4.5])
  C=[0.2;0.4;0.6;0.8]*[1 1 1]; colormap(C)
  [~,AA]=max(Q,[],2);
  patchplot(Xopt(:,2),Xopt(:,3),AA,1:4);
  xlabel('ponds')
  ylabel('adult ducks')
  title('Optimal Action')

