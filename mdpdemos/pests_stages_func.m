 clc
 disp('Stage models with P in memory and saved on disk')
 D=[0;5;20];  % damage costs
 C=10;        % spraying cost
 delta=0.975;  % discount factor
 T=inf;       % time horizon
 clear options
 X=[ ...
   1 1
   2 1
   3 1
   1 2
   2 2
   3 2];  
 % probabilities in growth phase 
 Gn=[ 0.65  0    0;
      0.25  0.55 0;
      0.10  0.45 1];
 Gt=[ 0.85  0    0;
      0.10  0.90 0;
      0.05  0.10 1];
 % probabilities in die-off phase
 Dn=[1  0.30 0.10;
     0  0.70 0.20;
     0  0    0.70];
 Dt=[1  0.80 0.50;
     0  0.20 0.30;
     0  0    0.20];
 Ix=getI(X,1);
 P1=[Gn Gt];
 pf{1}=saveget(P1);
 P2=[Dn Dt];
 pf{2}=saveget(P2);
 
   
 % set up model structure
 clear model
 model.name     = 'Stage model with P in memory';
 model.Ix       = Ix;
 model.ns       = 3;
 model.R        = {[zeros(3,1) zeros(3,1)-C];[-D -D-C]};
 model.discount = {1;delta};
 model.P        = {[Gn Gt],[Dn Dt]};
 model.T        = T;
 model.nstage   = 2;
 model.nrep = [1 1];
 % call solution procedure
 results=mdpsolve(model);
 v=results.v; x=results.Ixopt;
 disp(' ')
 disp('Staged Pest Control Problem')
 disp('State, Optimal Control and Value')
 fprintf('%2i  %2i  %2i %12.4f %12.4f\n',[(1:3)' X(x{1},2) X(x{2},2)  v{1} v{2}]')

 model.name = 'Stage model with P read from disk';
 model.P = pf;
 results=mdpsolve(model);
 v=results.v; x=results.Ixopt;
 disp(' ')
 disp('Staged Pest Control Problem')
 disp('State, Optimal Control and Value')
 fprintf('%2i  %2i  %2i %12.4f %12.4f\n',[(1:3)' X(x{1},2) X(x{2},2)  v{1} v{2}]')

disp(' ')
disp('Note: this demo saves information to the disk in the form of files with funny names')
disp('  and .MAT extensions. These can be erased. This uses the "saveget" procedure which')
disp('  does not have a method to clean up after itself')