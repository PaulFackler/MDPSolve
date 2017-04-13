% pest management problem when delta=1
 clc
 disp('Demonstrates the relative value and vanishing discount approaches to average reward problems')
 % create matrix of state/action combinations
 X=[ ...
   1 1
   2 1
   3 1
   1 2
   2 2
   3 2];
 D=[0;5;20];  % damage costs
 C=10;        % spraying cost
 % probabilities with no action
 P1=[0.65  0.15 0.05;
     0.25  0.40 0.20;
     0.10  0.45 0.75];
 % probabilities with spraying
 P2=[0.85 0.45 0.35;
     0.15 0.50 0.50;
     0.00 0.05 0.15];
 
 R=-[D D+C];
 
 % set up model structure
 clear model
 model.X        = X;
 model.discount = 1;
 model.R        = R;
 model.P        = [P1 P2];
 model.Ix       = [1 2 3 1 2 3];
 clear options
 options.algorithm = 'f';   % use policy iteration
 options.relval=1;
 options.tol=1e-8;
 options.getAopt=1;
 options.print=0;
 % call solution procedure
 results1=mdpsolve(model,options);
   
 disp('Simple Pest Control Problem')
 disp('State, Optimal Control and Value')
 disp('Relative value method')
 fprintf('%2i  %2i %10.6f\n',[results1.Xopt results1.v]')
 nu=results1.AR;
 disp(nu)
 
 options.vanish=0.9999;                 % use vanishing discount method
 options.relval=1;
 options.maxit=2000;
 results2=mdpsolve(model,options);
 fprintf('vanishing discount method with delta=%1.8f\n',options.vanish)
 nu2=results2.AR;
 fprintf('%2i  %2i %10.6f\n',[results2.Xopt results2.v]')
 disp(nu2)
 
 % check that the correct value is obtained by using the
 % long run distribution
 p=longrunP(results1.pstar);
 fprintf('Average reward computed using AR and P* %8.6f  %8.6f\n',nu,p'*R(results1.Ixopt))
 
 fprintf('maximum difference in value from the two algorithms: %1.6e\n',max(abs(results2.AR-nu)))
 
 
 options.vanish=0.9999;                 % use vanishing discount method
 options.relval=0;
 options.maxit=2000;
 options.algorithm='p';
 results3=mdpsolve(model,options);
 fprintf('vanishing discount method with delta=%1.8f\n',options.vanish)
 fprintf('%2i  %2i %10.6f\n',[results3.Xopt results3.v]')