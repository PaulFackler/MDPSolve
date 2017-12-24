% getprMC Uses Monte Carlo simulation to obtain reward and transition matrix
% USAGE
%   [P,R]=getprMC(D,options);
% INPUTS
%   D       : an influence diagram structure
%   options : structure variable with control options 
% OUTPUTS
%   P  : transition matrix
%   R  : reward function
%
% Options: see d2model (ptype is ignored as the transition matrix is
%            always computed)
%
% This procedure uses simulation to obtain estimates of the reward and 
% transition matrix. For each possible value of the states and actions
% it generates reps values of the future states and the utility function.
% R is equal to the mean value of the utility function conditional on the 
% current states and actions. P is computed using the mean of the interpolation
% weights of the future state values.

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

function [P,R]=getprMC(D,options)
reps            = 0;     % use Monte Carlo to construct reward and transition
                         % if reps>0 this is the # of replications to use
chunk           = 1000;  % chunk size for Monte Carlo approach
outvars         = [];    % output varariables
cleanup         = 0;     % used to handle extrapolation
sparseout       = false; % default is full
print           = 0;     % print level, 0: none, 1: moderate, 2: heavy
if nargin>=2 && ~isempty(options)
  if isfield(options,'reps'),        reps=options.reps;               end
  if isfield(options,'chunk'),       chunk=options.chunk;             end
  if isfield(options,'outvars'),     outvars=options.outvars;         end
  if isfield(options,'cleanup'),     cleanup=options.cleanup;         end
  if isfield(options,'sparseout'),   sparseout=options.sparseout;     end
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
values=D.values;

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
end

X=dvalues(D,[parameters actions states]);
nx=size(X{1},1);
s=values(outvars);
ns=prod(D.sizes(outvars));

evenspacing=zeros(1,length(s));
for i=1:length(s)
  ds=diff(s{i});
  if all(abs(ds-ds(1))<1e-14)
    evenspacing(i)=1;
  else
    evenspacing(i)=0;
  end
end

if ispc && ~sparseout 
  memstats=memory; 
  if memstats.MaxPossibleArrayBytes < ns*nx*8 
    sparseout=true; 
  end
end
if sparseout
  P=sparse(ns,nx);
else
  P=zeros(ns,nx);
end
R=zeros(nx,1);
j=1;
summat=kron(speye(chunk),ones(reps,1));
tstart=tic;
printed=false;
% loop over all state action combinations
while j<=nx
  if j+chunk-1>nx
    chunk=nx-j+1;
    summat=kron(speye(chunk),ones(reps,1));
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
  Sj=S(outvars);
  Pj=rectbas(Sj,s,evenspacing,cleanup)*summat;
  if issparse(P)
    P=add2sparse(P,Pj,j,[],true);
  else
    P(:,j:j+chunk-1)=Pj;
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
P=bsxfun(@rdivide,P,sum(P,1));  %P=P/reps;



% mdpsim  Monte Carlo simulation of an Markov process
% USAGE
%   ind = mdpsim(p,s)
% INPUTS
%   p   : nxm state transition matrix
%   s   : k by 1 vector of indices into the columns of p
%   z   : 1xk-vector of random uniform values
% OUTPUT
%   ind : k by 1 vector of indices to rows of p
function [ind,z] = mdpsim(p,s,z)
%n = size(p,1);
%u = ones(n,1);
k = length(s);
cp=cumsum(p,1); 
if nargin<3 || isempty(z)
  z = rand(1,k); 
else
  z=reshape(z,1,k);
end
%ind = 1+sum(r(u,:)>cp(:,s),1);
ind = 1+sum(bsxfun(@gt,z,cp(:,s)),1);


