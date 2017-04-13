% pest management problem - stage model
 clc
 disp('Pest management with growth and die-off stages')
 disp('Discounted case with no damages in growth stage')
 D1=[0;0;0];  % stage 1 damage costs
 D2=[0;5;20]; % stage 2 damage costs
 C=10;         % spraying cost
 delta=0.975;  % discount factor
 T=inf;        % time horizon
 clear options
 options.print=1;
 options.algorithm='f';
 options.tol=1e-12;
 options.maxit=2000;

 pests_stages_run
  
% get a case in which treatment in stage 1 is optimal
 disp('Discounted case with damages in growth stage')
 D1=D2;
 pests_stages_run

