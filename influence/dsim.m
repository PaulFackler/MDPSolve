% dsim Simulates variables in an Influence Diagram
% USAGE
%   [Y,z]=dsim(D,s0,T,A,pval,z,keepall);
% INPUTS
%   D       : an influence diagram structure
%   s0      : initial state values (reps x nd matrix)
%               if there are no state variables in D pass a scalar integer
%               representing the number of replications (reps)
%   T       : time horizon (positive integer)
%   A       : ns x da matrix representing the strategy or
%               a function handle of the form A(S1,S2,...,Snd)
%   pval    : parameter values (if any variables are parameter type an
%               assumed value must be specified)
%   z       : cell array of random values from previous call to dsim
%               for rvs with parents z contains that values of the underlying
%               random variables (uniform, normal or integer) that are used
%               to generate the rvs; for variables with no parents the values
%               of the rvs themselves are returned
%   keepall : d-element logical vector to indicate which variables are kept at 
%               every period; if keepall(i)=false only last period is returned
% OUTPUT
%   Y     : d-element cell array containing reps x T matrices, one
%             for each of the d variables in the diagram
%   z     : cell array of random values to use on subsequent calls to dsim
% 
% A(i,j) is the value of action j taken when the current state is state i. 
% The optimal strategy can be obtained using
%   Aopt=model.X(results.Ixopt,1:da);
% where model and results are the input and output of a call to mdpsolve
% and da is the number of action variables in the model

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2014-2017, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function [Y,z]=dsim(D,s0,T,A,pval,z,keepall)

d=length(D.names);
if ~any(ismember(D.types,'s'))
  reps=s0;
else
  reps=size(s0,1);
end
types=D.types;
if (nargin<5 || isempty(pval)) && any(ismember(types,'p'))
  error('parameter values must be specified when parameters are included')
end
if nargin<6, z=[]; end
% controls which variables are stored for all periods 
% (for others only the last period is kept)
if nargin<7 || all(keepall), keepall=true(1,d); end  
if isscalar(keepall), keepall=keepall & true(1,d); end

% determines whether z is actually kept
if nargout>1, keepz=true; else keepz=false; end  
cpds=D.cpds;
parents=getparents(D);
Y=cell(1,d);    % storage for all variables & time periods
St=cell(1,d);   % storage for all variables for current time period
ns=0;           % # of state variables
na=0;           % # of action variables
np=0;           % # of parameter variables
% get an index vector of state variables to use in determining
% the associated action. If A is numeric an index function
% is created to map values of S to rows of the S matrix (nearest neighbor)
% and the associated rows of A. 
svars=find(ismember(types,'s'));
if ~isempty(A)
  if isnumeric(A)
    statevals=cell(1,length(svars));
    for i=1:length(svars)
      statevals{i}=D.cpds{svars(i)}.values;
    end
    stateindexfunc=v2ifunc(statevals);
  end
end
fvars=find([types{:}]=='f');

%%%%% set up the simulation parameters for each variable
%%%%% initializing the states, parameters and actions
% if the ith variable is: 
%   an action match(i) is the associated column of A
%   a future state match(i) is the associated state variable
match=zeros(1,d);     
vartypes=zeros(1,d);  % 1 has no parents, 2 has parents, 3 has rvsimfunc
for i=1:d
  cpdi=cpds{i};
  if strcmp(cpdi.type,'d') && ~isempty(parents{i}) && isempty(cpdi.parameters)
    cpds{i}=getmatchfunc(cpdi,D.values(parents{i}));
  end
  if keepall(i), Y{i}=zeros(reps,T); end  % pre-allocate memory 
  switch types{i}
    case 's'  % state variable
      ns=ns+1;
      switch size(s0,1)
        case 1
          St{i}=s0(ones(reps,1),ns);
        case reps
          St{i}=s0(:,ns);
        otherwise
          error('S0 must have one or reps rows')
      end
      if keepall(i), Y{i}(:,1)=St{i}; end
    case {'a','d'}   % action variable
      na=na+1;
      match(i)=na;   % links ith variable to jth action
    case {'c','u','r','f','h'}
      if isfield(cpds{i},'valfunc')
        vartypes(i)=0;
      else
        if isempty(parents{i})
          if isfield(cpds{i},'simfunc'), vartypes(i)=1;
          else                           vartypes(i)=3;
          end
        else
          if isfield(cpds{i},'simfunc'), vartypes(i)=2;
          else                           vartypes(i)=4;
          end
        end    
      end
      if types{i}=='f'  % match future state to current state
        match(i)=find(ismember(D.names,D.names{i}(1:end-1)));
        St{i}=St{match(i)}; % set future state equal to associated current state
      end
    case 'p'   % parameter variable - these never change
      np=np+1;
      switch size(pval,1)
        case 1
          St{i}=pval(ones(reps,1),np);
        case reps
          St{i}=pval(:,np);
        otherwise
          error('pvals must have one or reps rows')
      end
    case 'n'
      % null variable passed forward
    otherwise
      error(['variable type ' types{i} ' not recognized'])
  end
end

% initialize cell array for random terms (z)
if isempty(z) 
  z=cell(1,d);    % storage for random noise terms for reuse
  for i=1:d
    switch types{i}
    case {'c','u','r','f','h'}
      z{i}=repmat({[]},1,T);
    end
  end
end

% loop over time periods
for t=1:T
  % set current states equal to previous future states
  if ~isempty(svars)
    St(match(fvars))=St(fvars);
  end
  % get the current actions as functions of current states
  if ~isempty(A)
    if isnumeric(A)
      %ind=gridmatch(St(stateind),statevals); % get the index values of the states
      aind=stateindexfunc(St{svars});
    else
      At=A([St{svars}]);
    end
  end
  for i=1:d
    switch types{i}
    case 's'                         % state variables
      % nothing to do
    case {'a','d'}                   % action variables
      if isnumeric(A)
        St{i}=A(aind,match(i));
      else
        St{i}=At(:,match(i));
      end
    case {'c','u','r','f','h'}       % everything else except parameters
      zit=z{i}{t};
      switch vartypes(i)
      case 0      % has evaluation function (valfunc)
        St{i}=cpds{i}.valfunc(St{parents{i}});
      case {1,2}  % has simulation function (simfunc)
        if isempty(zit), 
          switch cpds{i}.ztype
          case 'u'
            zit=rand(reps,1); 
          case 'n'
            zit=randn(reps,1); 
          case 'i'
            zit=randn(reps,1); 
          end
          % check for transformation function - used with discrete
          if isfield(cpds{i},'u2z');
            zit=cpds{i}.u2z(zit);
          end
        end
        if vartypes(i)==1               % no parents
          zit=cpds{i}.simfunc(zit);     % store the variable to avoid re-generating it
          St{i}=zit;  
        else                            % has parents
          St{i}=cpds{i}.simfunc(zit,St{parents{i}});
        end
      case 3                            % use rvgen - no parents
        zit=rvgen(reps,cpds{i},[],zit); % store the variable to avoid re-generating it
        St{i}=zit;
      case 4                            % use rvgen - has parents
        [St{i},zit]=rvgen(reps,cpds{i},St(parents{i}),zit);
      end
      if keepz, z{i}{t}=zit; end
    end
    if keepall(i)
      Y{i}(:,t)=St{i};
    end
  end
end
% only the last period value is returned when keepall(i) is false
for i=1:d
  if ~keepall(i), Y{i}=St{i}; end
end
