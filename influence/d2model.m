% d2model Defines an MDPSolve model structure from a diagram
% USAGE
%   [model,svars,xvars,DD]=d2model(D,options);
% INPUTS
%   D       : a diagram structure
%   options : structure variable with control options (described below)
% OUTPUTS
%   model : a model structure that can be passed to mdpsolve
%            (d field must be set before calling mdpsolve)
%   svars : cell array with the names of the state variables
%   xvars : cell array with the names of the state/action variables
%
% Options:
%   inc          : if unobserved states are present, this is the # of 
%                    subintervals in the belief state discretization
%   reps         : set to a positive integer to use Monte Carlo to construct
%                     reward and transition (equals # of MC replications)
%   chunk        : number of state/actions processed at a time using the
%                     Monte Carlo approach.
%   passforward  : 1 to pass functions forward [default]
%                  0 to replace functions with CPTs
%   cleanup      : cleanup variable passed to rectbas to handle extrapolation
%                    0) no adjustment 1) adjust at end 2) adjust for each shock
%   ptype        : 0 for conditional expectation function
%                  1 for transition probability matrix
%   order        :  the order of processing of the variables
%   orderalg     : algorithm to determine elimination order
%                     0 default hybrid
%                     1 greedy
%                     2 optimal
%   orderdisplay : 1 to print processing order information
%   forcefull    : 0 use sparse factors
%                  1 convert sparse factors to full
%   print          print level, 0: none, 1: moderate, 2: heavy
%   d            : discount factor for model
%   T            : time horizon for model
%   vterm        : terminal value function for model

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

function [model,svars,xvars,DD]=d2model(D,options)
% default option values
if nargin<2, options=[]; end

inc  = [];   % number of subintervals in the belief state discretization
print           = 0;    % print level, 0: none, 1: moderate, 2: heavy
if nargin>=2 && ~isempty(options)
  if isfield(options,'inc'),         inc=options.inc;                 end
  if isfield(options,'print'),       print=options.print;             end
end
obs=D.obs;
os=find(ismember(D.types,'s') & obs);  % unobserved states
us=find(ismember(D.types,'s') & ~obs); % unobserved states
op=find(ismember(D.types,'p') & obs);  % observed parameters
up=find(ismember(D.types,'p') & ~obs); % unobserved parameters
a=find(ismember(D.types,{'a','d'}));
if ~isempty(us) || ~isempty(up)
  y=find(ismember(D.types,'c') & obs);
  options.ptype=1; % if unobserved variables present get transition matrix
else
  y=[];
end


% no unobserved parameters or states 
if isempty(us) && isempty(up) %&& 0
  [model,DD,cost]=getmodel(D,[],options);
  svars=D.names{os};
  xvars=D.names{[a os]};
  if ~isempty(options)
    if isfield(options,'d'), model.d=options.d; end
    if isfield(options,'T'), model.T=options.T; end
    if isfield(options,'vterm'), model.vterm=options.vterm; end
  end
  return
end

pvars=[op up];
if isempty(pvars)
  np=0;
else
  pvals=dvalues(D,pvars);
  np=size(pvals{1},1);
end

x=[a os us];
z=zeros(1,length(x)-length(a));
k=0;
for i=1:length(os)
  k=k+1;
  z(k)=find(ismember(D.names,[D.names{os(i)} '+'])); 
end
for i=1:length(us)
  k=k+1;
  z(k)=find(ismember(D.names,[D.names{us(i)} '+']));
end
z=[z y];

DD=D;
parents=getparents(DD);
if np==0
  [model,DD,cost]=getmodel(DD,y,options);
else
  options.ptype=1;
  seed=randi(intmax,1,1);
  R=cell(1,np);
  P=cell(1,np);
  X=[];
  % get the models for each value of the parameters
  for i=1:np
    for j=1:length(pvars)
      DD.values{pvars(j)}=pvals{j}(i);
      DD.sizes(pvars(j))=1;
      DD.cpds{pvars(j)}=rvdef('d',1,pvals{j}(i));
    end
    for j=1:length(D.sizes)
      if any(ismember(parents{j},pvars)) && strcmp(DD.cpds{j}.type,'d');
        xx=dvalues(DD,parents{j});
        ind=gridmatch(xx,D.values(parents{j}));
        DD.cpds{j}.cpt=D.cpds{j}.cpt(:,ind);    
      end
    end
    rng(seed);
    [model,DD,cost]=getmodel(DD,y,options);
    P{i}=model.P;
    R{i}=model.R;
    X=[X;model.X];
    if np>0 && isfield(model,'Iexpand')
        P{i}=P{i}(:,model.Iexpand);
    end
    clear model
  end
end

x=[pvars x];

if isempty(up) && isempty(us)
  if np>0
    P=spblkdiag(P{:});
    R=[R{:}]; R=R(:);
    Ix=getI(X,[1:length(op) length(x)-length(os)+1:length(x)]);
    model=struct('P',P,'R',R(:),'X',X,'Ix',Ix);
  else
    model.Ix=getI(model.X,length(x)-length(os)+1:length(x));
  end
  xvars=D.names(x);
  svars=D.names(os);
else
  if isempty(inc)
    error('must specify inc option for POMDP models')
  end
  z=[pvars z];
  Z=dvalues(D,z,'m');

  if 1
    % determine the indices for the observed and unobserved states
    indox=[1:length(op) length(pvars)+length(a)+(1:length(os))];
    indux=[length(op)+(1:length(up)) length(pvars)+length(a)+length(os)+(1:length(us))];
    indoz=[find(ismember(z,op)) length(pvars)+1:length(pvars)+length(os)];
    induz=[find(ismember(z,up)) length(pvars)+length(os)+1:length(pvars)+length(os)+length(us)];
    if print>0
      disp(' ')
      disp('computing belief state model')
    end
    if np==0
      P=model.P;
      R=model.R;
      X=model.X;
    else
      P=spblkdiag(P{:});
      R=[R{:}]; R=R(:);
    end
    if exist('model','var') && isfield(model,'Iexpand'), P=P(:,model.Iexpand); end
    [Pb,Rb,Sb,Xb,Ixb]=xpomdp(inc,P,R,X,indox,indux,Z,indoz,induz);
  else
    if isempty(us) && isempty(y)  % try "standard" amdp if appropriate
      S=dvalues(D,[os us],'m');
      X=dvalues(D,[a os us],'m');
      Ix=getI(X,find(ismember([a os us],[os us])));
      if print>0
        disp('computing belief state model')
      end
      [b,Pb,Rb,Sb,Xb,Ixb]=amdp(inc,P,R,S,X,Ix);
    else                          % try "standard" POMDP if appropriate
      [A,AA]=adjacency(D);
      if ~(any(any(AA([os us],y))))  % y depends on future state
        if isempty(pvars)            % not set up to handle unknown parameters
          Q=conditional(D,y,[a os us]);
          [b,Pb,Rb]=pomdp(inc,[P{:}],Q,[R{:}],options);
          %Sb=rectgrid(dvalues(D,os),b);
          Xb=rectgrid(dvalues(D,[a os]),b);
          Ixb=getI(X,length(a)+1:length(a)+size(X,2));
        end
      end
    end
  end
  model=struct('P',Pb,'R',Rb,'X',Xb,'Ix',Ixb);
end
if ~isempty(options)
  if isfield(options,'d'), model.d=options.d; end
  if isfield(options,'T'), model.T=options.T; end
  if isfield(options,'vterm'), model.vterm=options.vterm; end
end

function [model,DD,cost]=getmodel(D,y,options)
reps            = 0;    % use Monte Carlo to construct reward and transition
                        % if reps>0 this is the # of replications to use
passforward     = 1;    % pass functions forward to children
cleanup         = 0;    % used to handle extrapolation
ptype           = 1;    % 0 to use condexp, 1 to use conditional
print           = 0;    % print level, 0: none, 1: moderate, 2: heavy
if nargin>=2 && ~isempty(options)
  if isfield(options,'reps'),        reps=options.reps;               end
  if isfield(options,'rep '),        reps=options.rep;                end
  if isfield(options,'passforward'), passforward=options.passforward; end
  if isfield(options,'cleanup'),     cleanup=options.cleanup;         end
  if isfield(options,'ptype'),       ptype=options.ptype;             end
  if isfield(options,'print'),       print=options.print;             end
end
  
states= find(ismember(D.types,'s'));
fstates=find(ismember(D.types,'f'));
actions=find(ismember(D.types,{'a','d'}));
rewards=find(ismember(D.types,{'r','u'}));
parameters= find(ismember(D.types,'p'));
if length(rewards)~=1
  error('A single reward (r) variable must be specified')
end
if length(states)+length(parameters)<1
  error('At least one state (s) variable or parameter variable (p) must be specified')
end
if length(fstates)<1
  %error('At least one future state (f) variable must be specified')
end
if length(actions)<1
  error('At least one action (d) variable must be specified')
end
if length(states)~=length(fstates)
  %error('There must be the same number of current (s) and future (f) state variables')
end

% ensure that current and future states are compatible and have the same order
fnames=D.names(fstates);
for i=1:length(fstates)
  fnames{i}(fnames{i}=='+')=[];
end
[junk,ii]=ismember(fnames,D.names(states)); %#ok<*ASGLU>
if any(ii==0)
  %error('Current (s) and future (f) state variables don''t match')
end
fnames=D.names(states);
for i=1:length(fstates)
  fnames{i}=[fnames{i} '+'];
end
[junk,fstates]=ismember(fnames,D.names);
xvars=[parameters(:); actions(:); states(:)]';

if reps==0
  if print>0
    disp(' ')
    disp('computing reward function')
    start=cputime;
  end
  uvar=find(ismember(D.types,{'u','r'}));
  options.getfunc=0;
  options.passforward=passforward;
  [R,~,~,Iexpand,DD]=condexp(D,uvar,xvars,options);
  if ~isempty(Iexpand)
    R=R(Iexpand);  
  end
  X=dvalues(D,xvars,'m');
  if any(isinf(R))
    feasible=R>-inf;
    R=R(feasible);
    X=X(feasible,:);
  else
    feasible=[];
  end
  
  model.R=R;
  model.X=X;
  model.Ix=getI(X,[1:length(parameters) length(actions)+1:length(xvars)]);
  
  if print>0
    fprintf('time taken to compute reward: %8.3f\n',cputime-start)
  end
  options.feasible=feasible;

  if ptype==1  % use transition matrix
    if print>0
      disp(' ')
      disp('computing transition matrix')
      start=cputime;
    end
    [model.P,V,cost,Iexpand] = conditional(DD,[fstates y],xvars,options);
    if ~isempty(Iexpand), 
      if length(Iexpand)<size(model.P,2)
        model.P=model.P(:,Iexpand);
      else
        model.Iexpand=Iexpand; 
      end
    end
    if print>0
      fprintf('time taken to compute transition matrix: %8.3f\n',cputime-start)
    end
  else  % use EV function
    if print>0
      disp(' ')
      disp('computing conditional expectation function')
      start=cputime;
    end
    options.getfunc=true;
    [model.P,V,cost]  = condexp(DD,fstates,xvars,options);
    model.EV = true;
    if print>0
      fprintf('time taken to compute conditional expectation function: %8.3f\n',cputime-start)
    end
  end
else % use Monte Carlo approach
  options.outvars=[fstates y];
  [model.P,R]=getprMC(D,options);
  X=dvalues(D,xvars,'m');
  if any(isinf(R))
    feasible=R>-inf;
    R=R(feasible);
    X=X(feasible,:);
    model.P=model.P(:,feasible);
  end
  model.R=R;
  model.X=X;
  model.Ix=getI(X,[1:length(parameters) length(actions)+1:length(xvars)]);
  DD=D;
  cost=0;
end

return
