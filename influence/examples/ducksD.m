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
clear variables
close all

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
            
Utility=@(H,AD) H.*min(1,AD/Ngoal); 
Utility=@(H,AD) H; 

% discretization for state variables
Nmin = .5;  Nmax = 18; Ninc = .5;
Pmin = .5;  Pmax = 8;  Pinc = .5;

% quadrature is used with the following number of nodes
np=11; % number of rain values 
nk=11; % number of harvest noise values
nn=11; % number of population noise values

N=(Nmin:Ninc:Nmax)';     % population in millions
P=(Pmin:Pinc:Pmax)';     % pond numbers in millions
HI=(1:4)';               % harvest intensity

% use RV structures for the noise terms
cptH = rvdef('n',[0;1]         ,nk); % harvest survival noise
cptN = rvdef('n',[0;sqrt(nvar)],nk); % population noise
cptR = rvdef('n',[0;sqrt(evar)],np); % rain noise
%%
D=add2diagram([],'ponds',             's',1,{}                                                                ,rvdef('v',[0,inf],P));
D=add2diagram(D, 'adult ducks',       's',1,{}                                                                ,rvdef('v',[0,inf],N));
D=add2diagram(D, 'harvest intensity', 'a',1,{}                                                                ,HI);
D=add2diagram(D, 'harvest noise',     'c',0,{}                                                                ,cptH);
D=add2diagram(D, 'population noise',  'c',0,{}                                                                ,cptN);
D=add2diagram(D, 'rain',              'c',1,{}                                                                ,cptR);
D=add2diagram(D, 'juvenile ducks',    'c',0,{'adult ducks','ponds'}                                           ,Repr);
D=add2diagram(D, 'harvest rate',      'c',0,{'harvest intensity','harvest noise'}                             ,HR);
D=add2diagram(D, 'ponds+',            'f',1,{'ponds','rain'}                                                  ,rvdef('f',Ptran,P));
D=add2diagram(D, 'adult ducks+',      'f',1,{'adult ducks','juvenile ducks','harvest rate','population noise'},rvdef('f',Nnext,N));
D=add2diagram(D, 'expected harvest',  'c',1,{'adult ducks','juvenile ducks','harvest intensity'}              ,rvdef('f',Harvest));
D=add2diagram(D, 'reward',            'r',1,{'expected harvest','adult ducks'}                                ,rvdef('f',Utility)); 

D.locs=[ ...
0.159 0.159 0.159 0.399 0.559 0.679 0.341 0.659 0.818 0.818 0.659 0.818;
0.798 0.599 0.299 0.898 0.898 0.898 0.699 0.449 0.798 0.599 0.299 0.080]';
D.attachments=[ ...
 2  1  3  4  1  6  2  7  8  5  2  7  3 11  2;
 7  7  8  8  9  9 10 10 10 10 11 11 11 12 12;
 6  4  5  4  5  5  5  5  6  4  4  4  5  5  4;
 2  8  1  1  1  8  1  1  2  8  8  8  1  7  1]';

figure(1); clf; set(gcf,'units','pixels','position',[680  450  1200 500]); drawdiagram(D)

%%
discretize=1; % 0 for simulation method

if discretize
  t=cputime;
  options=struct('ptype',1,'forcefull',0,'orderalg',0,'cleanup',0,'orderdisplay',1,'passforward',1,'print',1); 
  model=d2model(D,options);
  fprintf('time taken to set up model using d2model: %6.3f\n',cputime-t)
  if isnumeric(model.P)
    PP=model.P;
  else
    EV=model.P;
  end
else
  t=cputime;
  options=struct('reps',1000,'chunk',25,'cleanup',2,'print',1);
  model=d2model(D,options);
  fprintf('time taken to set up model using MC: %6.3f\n',cputime-t)
end
%%
model.d=0.98;
options.algorithm='i'; options.print=2;
results=mdpsolve(model,options);

%%
if isfield(results,'Ixopt') && ~isempty(results.Ixopt)
  Xopt=model.X(results.Ixopt,:);
  figure(2); clf
  set(gcf,'units','inches','position',[7  5  8 4.5])
  C=[0.2;0.4;0.6;0.8]*[1 1 1]; colormap(C)
  patchplot(Xopt(:,2),Xopt(:,3),Xopt(:,1),1:4);
  xlabel('ponds')
  ylabel('adult ducks')
  title('Optimal harvest regime strategy')
  legend('restricted','conservative','moderate','liberal','location','eastoutside')
else
  mdpreport(results)
end

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
