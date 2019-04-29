% Based on:
% U.S. Fish and Wildlife Service
% Adaptive Harvest Management, 2012 Hunting Season
% http://www.fws.gov/migratorybirds/NewReportsPublications/AHM/Year2012/AHMReport2012.pdf

% N:  population level
% P:  pond numbers
% HI: harvest intensity (1-4)
% HR: harvest rate
% HN: harvest noise
% H:  harvest
% AD: adult ducks
% JD: juvenile ducks
% K:  kill rate (HR/(1-c))

%close all
clear variables
close all
rng('default')

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

% pond transition function  (p. 34)
pP={'ponds','rain'};
fP = @(P,rain) max(0,2.2127 + 0.3420*P + rain);

% reproduction functions (p. 34)
pJD={'ponds','adult ducks','strongDD'};
fJD = @(P,AD,SR) ifthenelse(SR,1.1390 + 0.1376*P - 0.1131*AD, ...
                               0.7166 + 0.1083*P - 0.0373*AD);
pHR={'harvest intensity','harvest noise'};
fHR=@(HI,HN) HN.*sd(HI) + mk(HI);   % harvest rate
% survival functions
SurvA = @(K,sa) sa*max(0,(1-K));                            % additive survival
SurvC = @(K,sc) sc*(K<(1-sc)) + max(0,(1-K)).*(K>=(1-sc));  % compensatory survival
Surv=@(K,AS,s) ifthenelse(AS,SurvA(K,s),SurvC(K,s));

% duck transition function      
pAD = {'adult ducks','additiveH','juvenile ducks','harvest rate','population noise'};          
fAD = @(AD,AS,JD,h,PN) gams*AD.*(m*Surv(h,AS,s0ma)    + (1-m)*(Surv(afvr*h,AS,s0fa) + ...
                       gamr*JD.*(Surv(jfvr*h,AS,s0fa) + Surv(jmvr*h,AS,s0ma)*phi))).*exp(PN);
% function for the harvest level  
pEH = {'adult ducks','juvenile ducks','harvest intensity'};           
fEH = @(AD,JD,HI) gams*(m/phi+(1-m)*(afvr + (gamr*(jfvr+jmvr))*JD)).*mk(HI).*AD;

pU = {'expected harvest','adult ducks'}; 
fU = @(EH,AD) EH.*min(1,AD/Ngoal); 

% discretization for state variables
Pmin = .5;  Pmax = 8;  Pinc = .5;
Dmin = .5;  Dmax = 18; Dinc = .5;

P=(Pmin:Pinc:Pmax)';     % pond numbers in milliopns
AD=(Dmin:Dinc:Dmax)';    % population in millions
HI=(1:4)';               % harvest intensity

% quadrature is used with the following number of nodes
nk=5; % number of harvest noise values
nn=5; % number of population noise values
np=5; % number of rain values 

% use RV structures for the noise terms
HP =rvdef('n',[0;1]         ,nk); % harvest survival noise
PP =rvdef('n',[0;sqrt(nvar)],nn); % population noise
RP =rvdef('n',[0;sqrt(evar)],np); % rain noise

%%
p=cell(1,2); 
p{1}=g2P(fP,P,P,RP); 
g2=@(X,e)fAD(X(:,1),X(:,2),X(:,3),e(:,1),e(:,2));
p{2}=g2P(g2,AD,rectgrid(AD,AS,JD),{HP,PP});

%%
D=add2diagram([],'ponds',             's',1,{} ,P);
D=add2diagram(D, 'adult ducks',       's',1,{} ,AD);
D=add2diagram(D, 'harvest intensity', 'a',1,{} ,HI);
D=add2diagram(D, 'strongDD',          'p',0,{} ,[0;1]);
D=add2diagram(D, 'additiveH',         'p',0,{} ,[0;1]);
D=add2diagram(D, 'harvest noise',     'c',0,{} ,HP);
D=add2diagram(D, 'population noise',  'c',0,{} ,PP);
D=add2diagram(D, 'rain',              'c',0,{} ,RP);
D=add2diagram(D, 'juvenile ducks',    'c',0,pJD,fJD);
D=add2diagram(D, 'harvest rate',      'c',0,pHR,fHR);
D=add2diagram(D, 'ponds+',            'f',1,pP ,fP);
D=add2diagram(D, 'adult ducks+',      'f',1,pAD,fAD);
D=add2diagram(D, 'expected harvest',  'c',0,pEH,fEH);
D=add2diagram(D, 'reward',            'r',0,pU ,fU); 
D.locs=[ ...
0.120 0.120 0.120 0.249 0.380 0.415 0.599 0.731 0.434 0.636 0.878 0.878 0.620 0.878;
0.799 0.600 0.350 0.162 0.162 0.939 0.939 0.939 0.712 0.526 0.799 0.600 0.201 0.336]';
D.attachments=[ ...
 2 13  3  9  2  7 10  9  5  2  8  1  6  3  4  2  1;
14 14 13 13 13 12 12 12 12 12 11 11 10 10  9  9  9;
 4  5  5  4  4  4  5  5  6  5  4  5  4  5  6  6  4;
 1  1  1  8  1  8  1  1  2  1  8  1  8  1  2  1  1]';

figure(1); clf; set(gcf,'units','pixels','position',[100  100  1200 500]); drawdiagram(D)
%%
t=cputime;
options=struct('inc',5,'ptype',1,'orderalg',1,'cleanup',0,'reps',100,'print',0);
model=d2model(D,options);
fprintf('time taken to set up model: %6.3f\n',cputime-t)

%%
model.d=0.98;
options=struct('algorithm','i','print',2);
results=mdpsolve(model,options);

%%
if isfield(results,'Ixopt') && ~isempty(results.Ixopt)
  Xopt=model.X(results.Ixopt,:);
  figure(3-D.obs(4)); clf
  set(gcf,'units','inches','position',[3  3  8 6.5])
  C=[0.2;0.4;0.6;0.8]*[1 1 1]; colormap(C)
  pvals=dvalues(D,[4 5],'m');
  for i=1:4
    subplot(2,2,i)
    if D.obs(4)
      ind=Xopt(:,1)==pvals(i,1) & Xopt(:,2)==pvals(i,2);
      patchplot(Xopt(ind,4),Xopt(ind,5),Xopt(ind,3),1:4);
    else
      ind=Xopt(:,3+i)==1;
      patchplot(Xopt(ind,2),Xopt(ind,3),Xopt(ind,1),1:4);
    end
    xlabel('ponds')
    ylabel('adult ducks')
    title(['Optimal harvest regime strategy if B(' num2str(i) ') is true'])
  end
  %legend('restricted','conservative','moderate','liberal','location','eastoutside')
  legend('restricted','conservative','moderate','liberal','location','northeast')
else
  mdpreport(results)
end


