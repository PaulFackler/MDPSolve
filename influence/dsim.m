% dsim Simulates variables in an Influence Diagram
% USAGE
%   S=dsim(D,S0,reps,T,A,pval);
% INPUTS
%   D     : an influence diagram structure
%   S0    : initial state values (1 x ns vector or reps x ns matrix)
%   reps  : number of replications (number of paths)
%   T     : time horizon
%   A     : nsxda matrix representing the strategy
%   pval  : parameter values (if any variables are parameter type an
%             assumed value must be specified)
% OUTPUT
%   S    : d-element cell array containing reps x T+1 matrices, one
%            for each of the d variables in the diagram
% 
% A(i,j) is the value of action j taken when the current state is state i. 
% The optimal strategy can be obtained using
%   Aopt=model.X(results.Ixopt,1:da);
% where model and results are the input and output of a call to mdpsolve
% and da is the number of action variables in the model

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2014, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function S=dsim(D,S0,reps,T,A,pval)

d=length(D.names);
types=D.types;
if nargin<6 && any(ismember(types,'p'))
  error('parameter values must be specified when parameters are included')
end
cpds=D.cpds;
parents=getparents(D);
S=cell(1,d);
St=cell(1,d);
ns=0;
na=0;
np=0;
stateind=find(ismember(types,'s'));
statevals=cell(1,length(stateind));
for i=1:length(stateind)
  statevals{i}=D.cpds{stateind(i)}.values;
end
fvars=find([types{:}]=='f');

% set up the simulation parameters for each variable
% initializing the states, parameters and actions

% if the ith variable is 
%   an action match(i) is the associted column of A
%   a future state match(i) is the associated state variable
match=zeros(1,d);     
vartypes=zeros(1,d);  % 1 has no parents, 2 has parents
for i=1:d
  cpdi=cpds{i};
  if strcmp(cpdi.type,'d') && ~isempty(parents{i}) && isempty(cpdi.parameters)
    cpds{i}=getmatchfunc(cpdi,D.values(parents{i}));
  end
  S{i}=zeros(reps,T);
  switch types{i}
    case 's'
      ns=ns+1;
      switch size(S0,1)
        case 1
          St{i}=S0(ones(reps,1),ns);
        case reps
          St{i}=S0(:,ns);
        otherwise
          error('S0 must have one or reps rows')
      end
      S{i}(:,1)=St{i}; 
    case {'a','d'}
      na=na+1;
      match(i)=na;
    case {'c','u','r','f','h'}
      if isempty(parents{i})
        vartypes(i)=1;
      else
        vartypes(i)=2;
      end
      if types{i}=='f'
        match(i)=find(ismember(D.names,D.names{i}(1:end-1)));
      end
    case 'p'
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

% loop over time periods
for t=1:T
  if ~isempty(stateind)
    ind=gridmatch(St(stateind),statevals); % get the index values of the states
  end
  for i=1:d
    switch types{i}
    case 's'
      % nothing to do
    case {'a','d'}
      St{i}=A(ind,match(i));   
      S{i}(:,t)=St{i}; 
    case {'c','u','r','f','h'}
      switch vartypes(i)
      case 1
        St{i}=rvgen(reps,cpds{i});
      case 2
        St{i}=rvgen(reps,cpds{i},St(parents{i}));
      end
    end
    S{i}(:,t)=St{i};
    if types{i}=='f' && t<T
      S{match(i)}(:,t+1)=St{i};
    end
  end
  % change future states to current states
  if ~isempty(stateind)
    St(match(fvars))=St(fvars);
  end
end
return
