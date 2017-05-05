% dsim1step Uses Monte Carlo simulation to obtain reward and transition matrix
% USAGE
%   [S,R,X]=dsim1step(D,outvars,options);
% INPUTS
%   D       : an influence diagram structure
%   outvars : list of m variables to output
%   options : structure variable with control options 
% OUTPUTS
%   S  : m-element cell array of simulated values
%           Each element is reps x nx where nx is the number of 
%             combinations of state and action variables
%   R  : simulated reward function
%   X  : nx row matrix of parameter/action/state values
%
% Options: 
%   reps         : set to a positive integer to use Monte Carlo to construct
%                     reward and transition (equals # of MC replications)
%   chunk        : number of state/actions processed at a time using the
%                     Monte Carlo approach.
%
% This procedure generates reps simulated values of specified variables conditioned
% on all combinations of state and control values. 
% R is equal to the mean value of the utility function conditional on the 
% current states and actions. 

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

function [SS,R,X]=dsim1step(D,outvars,options)
if nargin<2, outvars=[]; end
reps            = 0;    % use Monte Carlo to construct reward and transition
                        % if reps>0 this is the # of replications to use
chunk           = 1000; % chunk size for Monte Carlo approach
print           = 1;    % print level, 0: none, 1: moderate, 2: heavy
if nargin>=2 && ~isempty(options)
  if isfield(options,'reps'),        reps=options.reps;               end
  if isfield(options,'chunk'),       chunk=options.chunk;             end
  if isfield(options,'print'),       print=options.print;             end
end
% get information on the variables
d=length(D.names);
states= find(ismember(D.types,'s'));
fstates=find(ismember(D.types,'f'));
actions=find(ismember(D.types,{'a','d'}));
parameters=find(ismember(D.types,'p'));
rewards=find(ismember(D.types,{'r','u'}));

parents=getparents(D);
cpds=D.cpds;

% determine variable types and initialize controls
match=zeros(1,d);
vartypes=zeros(1,d);
S=cell(1,d);
z=cell(1,d);  % storage for the underlying noise values
ns=0; 
na=0; 
np=0;
for i=1:d
  % get matchfunc for discrete rvs
  cpdi=cpds{i};
  if strcmp(cpdi.type,'d') && ~isempty(parents{i}) && isempty(cpdi.parameters)
    cpds{i}=getmatchfunc(cpdi,D.values(parents{i}));
  end
  switch D.types{i}
    case {'s'}
      ns=ns+1;
      match(i)=length(parameters)+length(actions)+ns;
    case {'a','d'}
      na=na+1;
      match(i)=length(parameters)+na;
    case 'p'
      np=np+1;
      match(i)=np;
    case {'c','u','r','f','h'}
      switch(cpdi.type)
        case 'f'
          vartypes(i)=2;
        otherwise
          if isempty(parents{i})
            vartypes(i)=1;
            S{i}=rvgen(reps,cpdi);  % pure noise - get once
          else
            z{i}=repmat(rvgen(reps,cpdi,[],'z'),chunk,1);
            vartypes(i)=3;
          end
      end
      if D.types{i}=='f'
        match(i)=find(ismember(D.names,D.names{i}(1:end-1)));
      end
    case 'n'
      % null variable passed forward
    otherwise
      error(['variable type ' D.types{i} ' not recognized'])
  end
end

% make the order of the future states the same
% as the order of the current states
[temp,ii]=sort(match(fstates)); fstates=fstates(ii); %#ok<ASGLU>
if isempty(outvars)
  outvars=fstates;
else
  if ~isnumeric(outvars)
    [junk,outvars]=ismember(outvars,D.names); %#ok<ASGLU>
  end
end

X=dvalues(D,[parameters actions states]);
nx=size(X{1},1);

R=zeros(nx,1);
j=1;
tstart=tic;
printed=false;
SS=cell(1,length(outvars));
for i=1:length(outvars)
  SS{i}=zeros(reps,nx);
end
% loop over all state action combinations
while j<=nx
  if j+chunk-1>nx
    chunk=nx-j+1;
  end
  for i=1:d
    switch D.types{i}
    case {'s','a','d','p'}
      S{i}=ones(reps,1)*X{match(i)}(j:j+chunk-1)';   
      S{i}=S{i}(:);
    case {'c','u','r','f'}
      switch vartypes(i)
      case 1  % pure noise case - already generated values
        if length(S{i})<reps*chunk
          S{i}=repmat(S{i},chunk,1);
        elseif length(S{i})>reps*chunk
          S{i}=S{i}((1:reps*chunk)');
        end
      case 2  % function value
        S{i}=double(cpds{i}.valfunc(S{parents{i}}));
      case 3  % defined by random variable structure
        if length(z{i})>reps*chunk 
          z{i}=z{i}(1:reps*chunk,1); 
        end
        S{i}=rvgen(reps*chunk,cpds{i},S(parents{i}),z{i});
      end
    end
  end
  for i=1:length(outvars)
    SS{i}(:,j:j+chunk-1)=reshape(S{outvars(i)},reps,chunk);
  end
  if ~isempty(rewards)
    R(j:j+chunk-1)=mean(reshape(S{rewards},reps,chunk),1); 
  end
  j=j+chunk;
  if print>=1 && toc(tstart)>5
    printed=true;
    if print==1
      fprintf('.')
    else
      fprintf('cases processed: %1i/%1i\n',j,nx)
    end
    tstart=tic;
  end
end
if printed==1, fprintf('\n'); end
if nargout>=3, X=[X{:}]; end

