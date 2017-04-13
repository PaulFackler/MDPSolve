% pest management problem
% multiple sites
function pests_catcount
 clc
 %close all
 disp('Multiple site pest infestation problem')
 pestscatcount(1)
 pestscatcount(0)
 
% if Pconstant=1, site probability does not depend on S
function pestscatcount(Pconstant)
 if nargin==0 || isempty(Pconstant), Pconstant=1; end
 prevent=1;   % if 1, allows treatment of uninfested sites
 
 N = 15;         % number of sites
 maxtreat = 5; % maximum number of treatments
 B = 60;        % annual spraying budget
 C=[5;15;20];   % treatment cost
 D=[0;20;50];    % damage costs
 delta=0.975;    % discount factor
 n=3;            % number of categories
 m=2;            % number of actions
 % probabilities with no action
 P1=[0.75  0.15 0   ;
     0.25  0.60 0.15;
     0     0.25 0.85];
 % probabilities with spraying
 P2=[0.90 0.65 0.35;
     0.10 0.30 0.55;
     0    0.05 0.10];
 P=[P1 P2];

 % X1-X3 are non treated, X4-X6 are treated
 X=simplexgrid(n*m,N,N,1,1,[N N N B./C']);
 ind=X*[zeros(n,1);C]<=B;  % identify the state/action combinations that do not exceed the budget
 if prevent==0
   ind=ind & X(:,4)==0;      % eliminates preventative treatment
 end
 ind=ind & sum(X(:,4:6),2)<=maxtreat+0.01;  % # of treatments must be less than maxtreat
 X=X(ind,:);
 [Ix,S]=getI(X*repmat(speye(n),m,1),1:n);  % X*repmat(speye(n),m,1) sums over all actions
 R=X*(-[D;D+C]);
 ns=size(S,1);
 
 if Pconstant
   disp('Constant site transition matrix')
   disp(P)
   tic
   Pcc=catcountP(N,3,6,P,X);
   fprintf('catcount time: %1.6f\n',toc)
 else
   b0u=log((1-P1(2,1))/P1(2,1));
   b0t=log((1-P2(2,1))/P2(2,1));
   b1u=zeros(1,6);
   b1t=zeros(1,6); 
   b1u([2 3])=[-0.5 -2.5];
   b1t([2 3])=[-0.25 -0.75];
   Pfunc=@(Xj) pests_catcountfunc(Xj,b0u,b0t,b1u,b1t,P1,P2,N);
   disp('Site transition matrix with no sites infested or treated')
   Pfunc([N 0 0 0 0 0]')
   disp('Site transition matrix with all sites severely infested and untreated')
   Pfunc([0 0 N 0 0 0]')
   tic
   Pcc=catcountP(N,3,6,Pfunc,X);
   fprintf('catcount time: %1.6f\n',toc)
 end
   
 % set up model structure
 clear model
 model.discount = delta;
 model.R        = R;
 model.P        = Pcc;
 model.Ix       = Ix;
 model.ns       = ns;
 model.colstoch = true;
 % call solution procedure
 tic
 options=struct('algorithm','p');
 results = mdpsolve(model,options);
 a=results.Ixopt; 
 toc
 
 numtreated=sum(X(a,[5 6]),2);
 mt=max(numtreated);
 
 if ns<20, shownum=true; else shownum=false; end
 if Pconstant
   figure(1); clf
   set(gcf,'name','Managing an infestation w/ transition probabilities constant')
 else
   figure(2); clf
   set(gcf,'name','Managing an infestation w/ transition probabilities dependent on severity of infestation')
 end
 C=linspace(0.9,0.3,mt+1)'*ones(1,3); colormap(C)  % uncomment for gray scale
 lnames=cell(1,mt+1);
 for i=0:mt
   lnames{i+1}=num2str(i);
 end
 subplot(1,3,3)
 h=zeros(1,3);
 h(1)=patchplot(S(:,2),S(:,3),numtreated,0:mt,shownum);
 axis square
 ylabel('# heavily infested')
 title('Total number of sites treated')
 for i=2:3
   subplot(1,3,+i-1)
   h(i+1)=patchplot(S(:,2),S(:,3),X(a,3+i),0:mt,shownum);
   axis square
   if i>1, xlabel('# moderately infested'); end
   if i==2, ylabel('# heavily infested'); end
   title(['Number of class ' num2str(i) ' sites treated'])
 end
 hh=legend(lnames,'orientation','horizontal');
 pos=get(hh,'position'); pos(2)=0.01; pos(1)=0.5-pos(3)/2; set(hh,'position',pos)
 set(gcf,'units','normalized','position',[.1 .25 .8 .5])
 
 