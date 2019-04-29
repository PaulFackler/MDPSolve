% Alabama deer model from Barry Grand
clear variables
%close all

%%%% Define variable names
% action
nHI='Harvest Intensity';
% current density states
nDD='Doe Density';
nDF='Fawn Density';
nDI='iBuck Density';
nDM='mBuck Density';
% harvest noise variables
nHND='Doe Harvest Noise';
nHNF='Fawn Harvest Noise';
nHNI='iBuck Harvest Noise';
nHNM='mBuck Harvest Noise';
% other mortality rates
nOMD='Doe Other Mortality';
nOMF='Fawn Other Mortality';
nOMI='iBuck Other Mortality';
nOMM='mBuck Other Mortality';
nFP ='Fawn Predation';
% harvest rates
nHRD='Doe Harvest Rate';
nHRF='Fawn Harvest Rate';
nHRI='iBuck Harvest Rate';
nHRM='mBuck Harvest Rate';
nH  ='Harvest';
% quality Metrics
nDEN='Total Density';
nAR ='Buck Age Ratio';
nSR ='Sex Ratio';
nDC ='Doe Condition';
nFEC='Fecundity';
nREC='Recruitment';
nQ  ='Quality';

hn = 7;         % number of noise values
dinc=2;
w = [1; 0.25; 0.75; 2]; % weights on harvest levels
w = [1; 1; 1; 1]; % weights on harvest levels

% harvest intensity actions
HI = (1:3)';
% density values
mDD = (1:dinc:25)';% Density Does					
mDF = (1:dinc:19)';% Density Fawns					
mDI = (1:dinc:11)';% Density Immature Bucks					
mDM = (1:dinc:25)';% Density Mature Bucks	

% harvest noise
[mHN,pHN]=qnwnormeven(hn,0,1);     % harvest survival noise
pHN=rvdef('n',[0;1],{mHN,pHN});
mk=[0.01 0.1 0.15]';              % mean harvest rates
sd=[0.001 0.01 0.02]';            % std. dev. of harvest rates
fHR=@(HN,HI) HN.*sd(HI) + mk(HI); % function for the harvest level  

fH=@(DD,DF,DI,DM,HRD,HRF,HRI,HRM) [DD.*HRD DF.*HRF DI.*HRI DM.*HRM]*w;
pH={nDD,nDF,nDI,nDM,nHRD,nHRF,nHRI,nHRM};

% other mortality rates
pOMD = rvdef('b',[1.8013;10.2075],hn);
pOMF = rvdef('b',[2.563;3.8444],hn);
pOMI = rvdef('b',[0.5118;2.9],hn);
pOMM = rvdef('b',[0.125;1.125],hn);
pFP  = rvdef('b',[4;6],hn);    % fawn predation
          
% functions for quality characteristics
fDEN = @(DD,DF,DI,DM) DD+DF+DI+DM;                                 % total deer density
fAR  = @(DI,DM)       DM./(DI+DM);                                 % age ratio of bucks
fSR  = @(DD, DI, DM)  (DI+DM)./(DD+DI+DM);                         % sex ratio of adults
fDC  = @(DD,DF,DI,DM) 0.3+(0.7./(1+exp(-((DD+DF+DI+DM)./45-1.25)./-0.2))); % doe condition (density dependent)
fFEC = @(DC)          (2.2./(1+exp(-(DC-0.7)./0.15)));             % fawns born
fREC = @(FP,FEC,DC)   ((1-FP).*FEC./(1+exp(-(DC-0.7)./0.09)));     % fawns recuited to fall population 
%fQ   = @(DEN, REC, AR, SR, DC) ...
%          DEN.*((0.2372.*REC./2.2)+(0.2308.*AR)+(0.2436.*abs(.4-SR)./.4)+(0.2885.*DC))./4; % overall quality
fQ   = @(DEN, REC, AR, SR, DC) ...
               ((0.2372.*REC./2.2)+(0.2308.*AR)+(0.2436.*abs(.4-SR)./.4)+(0.2885.*DC))./4; % overall quality
% parent names for quality characteristics        
pDEN = {nDD,nDF,nDI,nDM};
pAR  = {nDI,nDM};
pSR  = {nDD, nDI, nDM};
pDC  = {nDD,nDF,nDI,nDM};
pFEC = {nDC};
pREC = {nFP,nFEC,nDC};
pQ   = {nDEN, nREC, nAR, nSR, nDC};

% functions for future population states
fDD = @(DD, HRD, OMD, DF, HRF, OMF) DD.*(1-HRD).*(1-OMD)      + 0.5.*DF.*(1-HRF).*(1-OMF);
fDF = @(DD, HRD, OMD, REC)          DD.*(1-HRD).*(1-OMD).*REC;
fDI = @(DI, HRI, OMI, DF, HRF, OMF) 0.6.*DI.*(1-HRI).*(1-OMI) + 0.5.*DF.*(1-HRF).*(1-OMF);
fDM = @(DI, HRI, OMI, DM, HRM, OMM) 0.4.*DI.*(1-HRI).*(1-OMI) + DM.*(1-HRM).*(1-OMM);
% parent names for future population states
pDD = {nDD, nHRD, nOMD, nDF, nHRF, nOMF};
pDF = {nDD, nHRD, nOMD, nREC};
pDI = {nDI, nHRI, nOMI, nDF, nHRF, nOMF};
pDM = {nDI, nHRI, nOMI, nDM, nHRM, nOMM};

fU = @(Q,H) Q.*H; % make quality and weighted harvest the objective
pU = {nQ,nH};

%% create diagram
D=[];
% Initial population states
D=add2diagram(D,nDM,'s',1,{},mDM);
D=add2diagram(D,nDI,'s',1,{},mDI);
D=add2diagram(D,nDD,'s',1,{},mDD);
D=add2diagram(D,nDF,'s',1,{},mDF);
% action
D=add2diagram(D,nHI,'a',1,{},HI);

% harvest noise variables
D=add2diagram(D,nHNM,'c',0,{},pHN);
D=add2diagram(D,nHNI,'c',0,{},pHN);
D=add2diagram(D,nHND,'c',0,{},pHN);
D=add2diagram(D,nHNF,'c',0,{},pHN);

% other mortality rates
D=add2diagram(D,nOMM,'c',0,{},pOMM);
D=add2diagram(D,nOMI,'c',0,{},pOMI);
D=add2diagram(D,nOMD,'c',0,{},pOMD);
D=add2diagram(D,nOMF,'c',0,{},pOMF);
D=add2diagram(D,nFP ,'c',0,{},pFP );

% harvest rates
D=add2diagram(D,nHRM,'c',0,{nHNM,nHI},fHR);
D=add2diagram(D,nHRI,'c',0,{nHNI,nHI},fHR);  
D=add2diagram(D,nHRD,'c',0,{nHND,nHI},fHR);
D=add2diagram(D,nHRF,'c',0,{nHNF,nHI},fHR);

D=add2diagram(D,nH  ,'c',1,pH,fH); 

% Quality Metrics
D=add2diagram(D,nDEN,'c',1,pDEN,fDEN);
D=add2diagram(D,nAR ,'c',1,pAR ,fAR );
D=add2diagram(D,nSR ,'c',1,pSR ,fSR );
D=add2diagram(D,nDC ,'c',1,pDC ,fDC );
D=add2diagram(D,nFEC,'c',1,pFEC,fFEC);
D=add2diagram(D,nREC,'c',1,pREC,fREC);
D=add2diagram(D,nQ  ,'c',1,pQ  ,fQ  );

% Future population states 
D=add2diagram(D,[nDM '+'],'f',1,pDM,rvdef('f',fDM,mDM)); 
D=add2diagram(D,[nDI '+'],'f',1,pDI,rvdef('f',fDI,mDI)); 
D=add2diagram(D,[nDD '+'],'f',1,pDD,rvdef('f',fDD,mDD)); 
D=add2diagram(D,[nDF '+'],'f',1,pDF,rvdef('f',fDF,mDF)); 

% Reward
D=add2diagram(D,'Utility','r',0,pU,fU);

D.locs=[ ...
0.169 0.141 0.128 0.105 0.189 0.315 0.294 0.271 0.235 0.581 0.581 0.581 0.581 0.399 0.686 0.631 0.588 0.530 0.622 0.374 0.251 0.290 0.339 0.434 0.628 0.774 0.889 0.889 0.889 0.889 0.889;
0.669 0.600 0.529 0.460 0.740 0.948 0.898 0.850 0.800 0.421 0.370 0.321 0.273 0.418 0.948 0.898 0.850 0.800 0.708 0.324 0.147 0.208 0.265 0.061 0.061 0.061 0.669 0.600 0.529 0.460 0.303]';
D.attachments=[ ...
19 26 13 18  5 12 17  4 11 16  3 12 17  4 25 10 15  2 11 16  3 10 15  2 23 22 21 25 20 23 24 14 23  5  4  3  2  5  4  2  5  4  5  4  3  2 18 17 16 15  5  4  3  2  1  9  1  8  1  7  1  6;
31 31 30 30 30 30 30 30 29 29 29 29 29 29 28 28 28 28 27 27 27 27 27 27 26 26 26 26 26 25 25 25 24 23 23 23 23 22 22 22 21 21 20 20 20 20 19 19 19 19 19 19 19 19 18 18 17 17 16 16 15 15;
 4  6  5  5  5  5  5  5  5  5  5  5  5  5  6  5  5  5  5  5  5  5  5  5  5  5  5  5  5  4  5  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  6  6  6  6  5  5  5  5  5  5  5  5;
 1  1  2  8  1  2  8  1  2  8  1  2  8  1  2  2  8  1  2  8  1  2  8  1  8  8  8  1  8  1  1  8  1  8  8  8  8  8  8  8  8  8  8  8  8  8  8  8  8  8  1  1  1  1  1  1  1  1  1  1  1  1]';
%figure(1); set(1,'units','normalized','position',[.05 .1 .85 .8])
%drawdiagram(D,struct('fontsize',0.02))

%%
t=cputime;
options=struct('ptype',1,'forcefull',0,'orderalg',0,'orderdisplay',1,'print',1,'reps',1000); 
[model,svars,xvars]=d2model(D,options); 
fprintf('time taken to set up model using d2model: %6.3f\n',cputime-t)
if isnumeric(model.P)
  model.P=sparse(model.P);
  P=model.P;
else
  EV=model.P;
end
return


model.d=0.98;
options=struct('algorithm','i','print',2);
results=mdpsolve(model,options);
Aopt=model.X(results.Ixopt,1);

%% simulate the model with low intensity and optimal intensity
ns=size(dvalues(D,find(ismember(D.types,'s')),'m'),1);
rep=1000;
T=100;
S0=[17 9 5 15];
S=dsim(D,S0,rep,T,ones(ns,1));
m=[mean(S{2}(:,end),1) mean(S{3}(:,end),1) mean(S{4}(:,end),1) mean(S{5}(:,end),1)];
disp('Long-run mean population levels with low harvest action')
disp(m)

S=dsim(D,S0,rep,T,Aopt);
m=[mean(S{2}(:,end),1) mean(S{3}(:,end),1) mean(S{4}(:,end),1) mean(S{5}(:,end),1)];
disp('Long-run mean population levels with optimal harvest action')
disp(m)

%%
Xopt=model.X(results.Ixopt,:);
mu=Xopt(match(m,Xopt(:,2:5)),2:5);
figure(3); clf
set(gcf,'units','pixels','position',[300 75   1600  900])
C=[0.2; 0.45; 0.7]*[1 1 1]; colormap(C)
nDF=length(mDF);
nDI=length(mDI);
for i=1:nDI
  for j=1:nDF
    subplot(nDI,nDF,j+(nDI-i)*nDF)
    ind=Xopt(:,3)==mDF(j) & Xopt(:,4)==mDI(i);
    patchplot(Xopt(ind,2),Xopt(ind,5),Xopt(ind,1),[1 2 3]);
    if i==nDI, title(['Fawns=' num2str(mDF(j))]); end
    if j==1, ylabel(['iBucks=' num2str(mDI(i))]); end
    if mu(2)==mDF(j) && mu(3)==mDI(i) 
      hold on; plot(mu(1),mu(4),'ws'); hold off
    end
  end
end
h=legend('restrictive','moderate','liberal');
pos=get(h,'position'); pos(1)=1-pos(3)*1.05; pos(2)=0.5-pos(4)/2; set(h,'position',pos);


%% approximate the optmal strategy by fitting a multinomial logit model
mlresults = multilogit(Xopt(:,1)-1,[ones(size(Xopt,1),1) Xopt(:,2:5) Xopt(:,2:5).^3]);
[vv,Afit]=max(mlresults.yfit,[],2);
fprintf('concordence between optimal and fit: %6.3f\n', sum(Afit==Xopt(:,1))/size(Xopt,1))

figure(4); clf
set(gcf,'units','pixels','position',[300 75   1600  900])
C=[0.2; 0.45; 0.7]*[1 1 1]; colormap(C)
nDF=length(mDF);
nDI=length(mDI);
for i=1:nDI
  for j=1:nDF
    subplot(nDI,nDF,j+(nDI-i)*nDF)
    ind=Xopt(:,3)==mDF(j) & Xopt(:,4)==mDI(i);
    patchplot(Xopt(ind,2),Xopt(ind,5),Afit(ind),[1 2 3]);
    if i==nDI, title(['Fawns=' num2str(mDF(j))]); end
    if j==1, ylabel(['iBucks=' num2str(mDI(i))]); end
    ind=ind&Afit~=Aopt(:,1);
    hold on; plot(Xopt(ind,2),Xopt(ind,5),'wx'); hold off
    if mu(2)==mDF(j) && mu(3)==mDI(i) 
      hold on; plot(mu(1),mu(4),'ws'); hold off
    end
  end
end
h=legend('restrictive','moderate','liberal');
pos=get(h,'position'); pos(1)=1-pos(3)*1.05; pos(2)=0.5-pos(4)/2; set(h,'position',pos);
