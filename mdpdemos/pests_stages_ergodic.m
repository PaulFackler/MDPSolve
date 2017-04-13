% pest management problem - stage model
 clc
 disp('Pest management with growth and die-off stages - average reward case')
 D1=[0;0;0];  % stage 1 damage costs
 D2=[0;5;20]; % stage 2 damage costs
 % uncomment next line to get a case in which treatment in stage 1 is optimal
 D1=D2;
 C=10;        % spraying cost
 delta=1;     % discount factor
 T=inf;
 
 clear options
 options.print=0;
 options.relval=true;
 options.algorithm='p';
  
 pests_stages_run