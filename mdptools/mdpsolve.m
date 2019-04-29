% mdpsolve  Solves discrete-state/action dynamic programs
% USAGE
%   results = mdpsolve(model,options);
% INPUTS
%   model   : model structure variable (see below)
%   options : options structure variable (see below)
% OUTPUTS
%   results : structure variable with problem results
%               v           value function
%               Ixopt       index of optimal rows of X
%               Aopt        optimal actions (if requested with restrictions)
%               pstar       optimal transition probability matrix (only available if T==inf)
%               algorithm   'p', 'f' or 'b'
%               time        run time
%               iter        number of iterations (for T=inf only)
%               change      maximal change on the last iteration (for T=inf only)
%               numnochange number of interations since the last change in policy
%               errors      cell array containing error information
%               warnings    cell array containing warning information
% Error and warning information is printed if options.print>=1
% It can also be viewed using mdpreport(results)
%
% For non-stage models v and Ixopt are either n x 1 or 
%   are n x T if T<inf and keepall=1.
%   The optimal strategy is found using X(Ixopt,:) or X(Ixopt(:,t),:).
% For stage models v and Ixopt are cell arrays with nstage elements
%   with the ith element n x nrep(i) x 1 or n x nrep(i) x T.
%   The optimal strategy is found using X(Ixopt{i}(:,j,:) or X(Ixopt{i}(:,j,t),:)
%     for stage i, replication j.
%
% The possible fields of the model structure variable are:
%   d or discount  : a scalar on (0,1]
%   R or reward    : nx-vector of reward values
%   transfunc      : for deterministic problems 
%                      an nxx1 vector of values on {1,...,ns} (not currently implemented)
%   P or transprob : There are several ways to define the transition matrices
%                      an ns x nx array (if colstoch=1)
%                      an nx x ns array (if colstoch=0)
%   T              : number of time periods 
%                     (omit or set to inf for infinite horizon problems)
%   vterm          : terminal (time T+1) value function
%   Ix             : nx-vector indicatign which state value is associated with
%                      each state/action combination
%   Iexpand        : nx-vector of indices specifying the column of P associated
%                      with each state/action combination (significant
%                      memory and speed efficiencies possible using this option)
%   EV             : 0/1 indicating that P is a function that accepts ns-vector V 
%                      and returns nx vector E[V+|X]
%   colstoch       : 0/1 indicating that P is column stochastic (rows are future
%                      values). This is generally not needed unless the Iexpand
%                      feature is used because the number of states is less than
%                      the number state/action combinations
%   X              : nx-row matrix of all state/action values
%   svars          : which columns of X are state variables
% R, d and P must be specified. If T is omitted, an infinite horizon problem is assumed. 
% If vterm is omitted, the terminal value if assumed to be 0. If Ix is
% omitted R should be an ns row matrix with nx elements (implying that there
% are the same number of actions for each state) or X and svars should be specified.
%
% The possible fields of the options structure variable are:
%       FOR ALL HORIZON PROBLEMS
%   print       : 0 no output to the screen
%                 1 display summary report
%                 2 display information at each iteration and summary report
%   checks      : 0/1 checks perfoms checks on input data
%       FOR FINITE HORIZON PROBLEMS
%   keepall     : keep values and actions for every iteration
%       FOR INFINITE HORIZON PROBLEMS
%   algorithm   : 'p' for policy iteration, 'f' for function iteration
%                   for infinite horizon problems (default: 'policy')
%   modpol      : non-negative integer equal to the maximum number of times to
%                   run modified policy iterations (for algorithm='f' only)
%                   default: 100
%   relval      : controls the use of the relative value algorithm
%                  0: uses ordinary policy or function iteration
%                  relval = k in {1,2,...,ns} uses the relative value algorithm 
%                  with the value relative to state k, i.e., V(k)=0 is set to 0.
%                  In this case the value of V(k) is returned as results.AR
%   vanish      : controls the use of the vanishing discount approach. Set to 0
%                 to not use this, set to a number close to but less than 1 to
%                 obtain a discounted approximation to the non-discounted case.
%   maxit       : maximum number of iterations
%   tol         : convergence tolerance
%   nochangelim : stop if action does not change in nochangelim iterations
%   v           : starting value vector
%
% Other options are available for specifying a model. See user documentation.

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2013, Paul L. Fackler (paul_fackler@ncsu.edu)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without  
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%        this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice, 
%        this list of conditions and the following disclaimer in the 
%        documentation and/or other materials provided with the distribution.
%    * Neither the name of the North Carolina State University nor of Paul L. 
%        Fackler may be used to endorse or promote products derived from this 
%        software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% For more information, see the Open Source Initiative OSI site:
%   http://www.opensource.org/licenses/bsd-license.php

function results = mdpsolve(model,options)
% EXTRACT INFORMATION FROM MODEL STRUCTURE
 tic;
% Default options
  print       = 1;     % print interation information to screen
  checks      = 1;     % performs consistency checks (set to 0 it avoid these checks)
  keepall     = 0;     % keep values and actions for every iteration (T<inf only)
  algorithm   = 'p';   % algorithm (p: policy iteration, f: function iteration)
  modpol      = 100;   % number of times to run modified policy iterations (for algorithm='f' only)
  relval      = false; % state used with the relative value approach
  vanish      = false; % set to a number close to but less than 1 to obtain the
                       % vanishing discount method (when delta=1)
  maxit       = [];    % set below                      
  tol         = 1e-8;  % convergence tolerance (for algorithm='f' only)
  nochangelim = inf;   % stop if action does change in nochangelim iterations
                       % (for algorithm='f' only)
  v           = [];    % starting value vector (T<inf only)
  debug       = 0;
  % set default maximu number of iterations

  if nargin>=2 && ~isempty(options)
    if isfield(options,'print'),       print=options.print;             end
    if isfield(options,'checks'),      checks=options.checks;           end
    if isfield(options,'keepall'),     keepall=options.keepall;         end
    if isfield(options,'algorithm'),   algorithm=options.algorithm;     end
    if isfield(options,'modpol'),      modpol=options.modpol;           end
    if isfield(options,'maxit'),       maxit=options.maxit;             end
    if isfield(options,'relval'),      relval=options.relval;           end
    if isfield(options,'vanish'),      vanish=options.vanish;           end
    if isfield(options,'tol'),         tol=options.tol;                 end
    if isfield(options,'nochangelim'), nochangelim=options.nochangelim; end
    if isfield(options,'v'),           v=options.v;                     end
    if isfield(options,'debug'),       debug=options.debug;             end
  end
  
 [errors,warnings,R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,T,nstage,nrep,Xindexed,expandP] ...
    = mdp_unpack(model,debug);

  % if errors occurred in unpacking the model, exit now
  if ~isempty(errors)
   results.errors=errors; 
   if ~isempty(warnings), results.warnings=warnings; end
   if print>0,  mdpreport(results); end
   return
 end
 warn0={};
    
  algorithm = lower(algorithm(1)); 
  if algorithm=='n', algorithm='p'; end  % change newton to policy
  
  if isempty(maxit)
    if algorithm=='p',   maxit = ceil(20*log(double(ns(1))));         
    else                 maxit = ceil(250*log(double(ns(1))));
    end
  end
    
  % determine appropriate algorithm (functionit, policyit or backit)
  % and get terminial value or starting value
  if T<inf
    if isfield(model,'vterm'),  v = model.vterm(:);
    else                        v=zeros(ns(1),1);
    end
  else
    if ~(algorithm =='p' || algorithm=='f' || algorithm=='i')
      results.errors={31};  % incorrect algorithm choice
      return
    end
    if algorithm=='p' %&& any(EV)
      %warn0{1,end+1}={51}; % Policy iteration not implemented with EV option
      if isfield(model,'d'), discount=model.d;
      else                   discount=model.discount;
      end
      if length(discount)==1 && discount<1
        algorithm='i';
      else
        algorithm='f';
      end
    end
    if isempty(v),  v=zeros(ns(1),1); 
    else            v=v(:);    
    end
  end
  % check size of initial value function vector
  if length(v)~=ns(1)
    if T<inf
      results.errors={32};
    else
      results.errors={33};
    end
    return
  end
  
  % call appropriate solver
  if nstage==1 
    [errors,warnings,R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP] =  ...
        mdp_getparams(1,checks,R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP,debug);
    if ~isempty(errors)>0
      results.errors=errors; results.warnings=[warn0 warnings];
    else
      if T==inf  % infinite horizon, non-stage model
        results = mdpsolve_Inf(R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP, ...
          v,algorithm,modpol,relval,vanish,maxit,tol,nochangelim,print,[]);
      else       % finite horizon, non-stage model
        results = mdpsolve_Fin( ...
           R, P, d, ns, nx, Ix, Iexpand, colstoch, EV, Xindexed, expandP, T, v, keepall, print);
      end
    end
  else
    if T==inf  % infinite horizon, stage model
      results = mdpsolve_Inf_stage(R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP, ...
         nstage,nrep,v,algorithm,modpol,relval,vanish,maxit,tol,nochangelim,print,checks);
    else       % finite horizon, stage model
      results = mdpsolve_Fin_stage(nstage,nrep, ...
         R, P, d, ns, nx, Ix, Iexpand, colstoch, EV, Xindexed,expandP, T, v, keepall, print,checks,debug);
    end 
  end
  % combine the warning generated here, from mdp_getparams and from solver
  if ~isfield(results,'warnings')
    results.warnings=[warn0 warnings];
  else
    results.warnings=[warn0 warnings results.warnings];
  end
  if ~isfield(results,'errors')
    results.errors={};
  end
  % for non-staged models get optimal variable values
  results.Xopt=[];
  if nstage==1 && isfield(model,'X')
    try
      results.Xopt=getA(model.X,results.Ixopt);
    catch ME
      results.warnings{end+1}={97,ME.message};
    end
  end
  % reorder the fields
  if isfield(model,'name')
    results.name=model.name;
  end
  if T<inf
    results.algorithm='b';
  end
  fnames={'name','v', 'AR','Ixopt', 'Xopt', 'pstar','algorithm','time','iter','stage',...
          'MPI','change','numnochange','errors','warnings'};
  for i=1:length(fnames)
    if ~isfield(results,fnames{i})
      results.(fnames{i})=[];
    end
  end
  results=orderfields(results,fnames);
  results.time=toc;
  if print>0
    mdpreport(results)
  end
  
  