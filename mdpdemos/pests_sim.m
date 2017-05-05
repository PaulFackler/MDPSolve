 clc

% pest management problem
% demonstrates simulation procedure

 D=[0;5;20];  % damage costs
 C=10;        % spraying cost
 delta=0.95;  % discount factor
% probabilities with no action
 P1=[0.65  0.15 0.05;
     0.25  0.40 0.20;
     0.10  0.45 0.75];
% probabilities with spraying
 P2=[0.85 0.45 0.35;
     0.15 0.50 0.50;
     0    0.05 0.15];

 R=-[D D+C];
 P=[P1 P2];

% set up model structure
 clear model
 model.discount = delta;
 model.R        = R;
 model.P        = P;
 clear options
 options.algorithm = 'p';
% call solution procedure
 results=mdpsolve(model,options);
 v=results.v; x=results.Ixopt; pstar=results.pstar;

 T=25;
 St=mdpsim((1:3)',P,x,T,true);

 figure(1); clf
 h=plot(0:T,St,'.-');
 title('Simulated Time Paths of State for Pest Problem')
 xlabel('time')
 ylabel('S')
 ylim([0.75 3.25])
 set(gca,'ytick',(1:3));

 T=5000;
 St=mdpsim((1:3)',P,x,T,true);
 pp=zeros(3,1);
 for i=1:3
   pp(i)=sum(St(:)==i);
 end
 lrp=longrunP(P(:,x));
 disp(['Exact and Estimated Long-run Probabilities (T=' num2str(T) ')'])
 disp([lrp pp/(3*T)])

