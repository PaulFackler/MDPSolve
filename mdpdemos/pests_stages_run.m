% pest management problem - stage model
% should be called after setting parameter values and mdpsolve options
% see, e.g., pests_stages
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
   
 % set up model structure
 clear model1
 Ix=getI(X,1);
 model1.name     = 'Model 1';
 model1.Ix       = Ix;
 model1.ns       = 3;
 model1.R        = {[-D1 -D1-C],[-D2 -D2-C]};
 model1.discount = {1;delta};
 model1.discount = sqrt(delta);
 model1.P        = {[Gn Gt],[Dn Dt]};
 model1.T        = T;
 model1.nstage   = 2;
 model1.nrep     = [1 1];
 model1.colstoch=true;
 % call solution procedure
 results1=mdpsolve(model1,options);
 v=results1.v; x=results1.Ixopt;
 disp(' ')
 disp('Staged Pest Control Problem - using MDPSOLVE stage feature')
 disp('State, Optimal Control and Value')
 fprintf('%2i  %2i  %2i %12.4f %12.4f\n',[(1:3)' X(x{1},2) X(x{2},2)  v{1} v{2}]')
 if delta==1
   fprintf('average reward: %12.4f\n',results1.AR); 
 end

 % Action in Stage 2 only
 % set up model structure
 XX        = {X(1:3,:) X};
 Ix1=getI(XX{1},1);
 Ix2=getI(XX{2},1);
 model2=model1;
 model2.name     = 'Model 2';
 model2.Ix       = {Ix1,Ix2};
 model2.R        = {-D1; [-D2 -D2-C]};
 model2.P        = {Gn,[Dn Dt]};
 % call solution procedure
 results2=mdpsolve(model2,options);
 v=results2.v; x=results2.Ixopt;
 disp('Staged Pest Control Problem - no treatment allowed at first stage')
 disp('State, Optimal Control and Value')
 fprintf('%2i  %2i  %2i %12.4f %12.4f\n',[(1:3)' XX{1}(x{1},2) XX{2}(x{2},2)  v{1} v{2}]')
 if delta==1
   fprintf('average reward: %12.4f\n',results2.AR); 
 end
 
 % Non-staged model
 % set up model structure
 Ix=getI(X,1);
 clear model3
 model3.name     = 'Model 3';
 model3.Ix       = Ix;
 model3.R        = -[D2+sqrt(delta)*Dn'*D1 D2+sqrt(delta)*Dt'*D1+C];
 model3.discount = delta;
 model3.P        = [Gn*Dn Gn*Dt];
  % call solution procedure
 results3=mdpsolve(model3,options);
 v=results3.v; x=results3.Ixopt;
 disp('Non-staged Pest Control Problem - action taken after growth phase only')
 disp('State, Optimal Control and Value')
 fprintf('%2i  %2i  %12.4f\n',[(1:3)' X(x,2) v]')
 if delta==1
   fprintf('average reward: %12.4f\n',results3.AR); 
 end
 
 X=rectgrid([1;2],X);
 Ix=getI(X,[1 2]);
 clear model4
 model4.name     = 'Model 4';
 model4.R        = [-D1; -D1-C; -D2; -D2-C];
 model4.discount = sqrt(delta);
 model4.P        = [zeros(3,3) zeros(3,3) Dn Dt;Gn Gt zeros(3,3) zeros(3,3)];
 model4.Ix=Ix;
 model4.ns=6;
 model4.colstoch=true;
 options.nochangelim=inf;
 results4=mdpsolve(model4,options);
 v=results4.v; x=results4.Ixopt;
 disp('Non-staged Pest Control Problem - solve stages simultaneously')
    disp('stage state  action     value')
 fprintf(' %2i    %2i      %2i     %12.4f\n',[X(x,:) v]')
 if isfield(options,'relval') && options.relval>=1 && delta==1
   disp('The relative value algorithm may give incorrect results in this case due to periodicity')
 end
 if delta==1
   fprintf('average reward: %12.4f\n',2*results4.AR); 
 end
 