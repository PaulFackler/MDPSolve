% condexp Creates a function to compute conditional expectations
% This creates a function to compute Ef=E[f(children)|parents],
%   a function that takes an nc vector and returns an np vector.
% It must be the case that the parents are not descendants of the children.
% USAGE
%   [Ef,Iexpand,V,cost,DD]=condexp(D,clist,plist,options);
% INPUTS
%   D      : an influence diagram (created using add2diagram)
%   clist  : the variables (children) over which expectations are computed
%              passed as a vector of variable numbers or a cell array of
%              diagram names (strings)
%   plist  : the conditioning variables (parents) passed as a vector of variable 
%              numbers or a cell array of diagram names (strings)
%   options : structure variable with control options (described below)
% OUTPUTS
%   Ef       : a function handle for a function of the form Ef(x)
%                where x is a vector of length nc=size(dvalues(D,clist,'m'),1)
%   V        : (m-1)x7 cell array containing information on the order
%                factors are processed and the variables involved
%                Use with orddisp function.
%   cost     : total processing cost (number of multiply operations)
%   Iexpand  : If requested Iexpand is an index vector of length
%                np=size(dvalues(D,plist,'m'),1) that can be used to expand
%                the output of Ef if some of the variables in plist are
%                not ancestors of the clist variables.
%                Thus 
%                  efx=Ef(x); efx=efx(Iexpand);
%                has the proper size. If not requested Iexpand is applied
%                in the Ef function.
%
% Options:
%   passforward  : 1 to pass functions forward [default]
%                  0 to replace functions with CPTs
%   cleanup      : cleanup variable passed to rectbas to handle extrapolation
%                    0) no adjustment 1) adjust at end 2) adjust for each shock
%   order        :  the order of processing of the variables
%   orderalg     : algorithm to determine elimination order
%                     0 default hybrid
%                     1 greedy
%                     2 optimal
%   orderdisplay : 1 to print processing order information
%   forcefull    : 0 use sparse factors
%                  1 convert sparse factors to full
%                  2 convert sparse factors used in EV function to full 
%   feasible     : logical vector for the state/action values that are feasible
%                    (R>-inf)

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

function [Ef,V,cost,Iexpand,DD]=condexp(D,clist,plist,options)
% default option values
if nargin<4, options=[]; end
passforward     = 1;   % pass functions forward to children
cleanup         = 0;   % used to handle extrapolation
forcefull       = 0;   % 1 to force sumproduct to use full factors
orderdisplay    = 0;   % displays sum-product order information
getfunc         = 0;   % 1 to get a function handle
if nargin>=4 && ~isempty(options)
  if isfield(options,'passforward'),  passforward=options.passforward;   end
  if isfield(options,'cleanup'),      cleanup=options.cleanup;           end
  if isfield(options,'forcefull'),    forcefull=options.forcefull;       end
  if isfield(options,'orderdisplay'), orderdisplay=options.orderdisplay; end
  if isfield(options,'getfunc'),      getfunc=options.getfunc;           end
end

% future states are default children  
if nargin<2 || isempty(clist)
  clist=find(ismember(D.types,{'f'}));
end
% current states, actions and parameters are default parents
if nargin<3 || isempty(plist)
  plist=find(ismember(D.types,{'s','d','a','p'}));
end
if ~isnumeric(clist), [junk,clist]=ismember(clist,D.names); end %#ok<*ASGLU>
if ~isnumeric(plist), [junk,plist]=ismember(plist,D.names); end

if length(clist)>1 && ~getfunc
  error('clist must be a singleton if getfunc option is false')
end

DD=D;
DD.parents=getparents(DD); % convert names to numbers
for i=1:length(clist)
  var=clist(i);
  cpd=DD.cpds{var};
  if strcmp(cpd.type,'f') 
    if ~isfield(cpd,'expvalues')
      [X,pvars,DD]=dvalues(DD,DD.parents{var},options);
      cpd.expvalues=double(cpd.valfunc(X{:}));
      cpd.expparents=pvars;
      DD.cpds{var}=cpd;
    end
    if getfunc 
      if ~isfield(cpd,'cpt') 
        DD=f2cpt(DD,var,cleanup,passforward);
      end
    else
      DD.values{var}=cpd.expvalues;
      DD.parents{var}=cpd.expparents;
      DD.sizes(var)=length(cpd.expvalues);
    end
  end
end



[A,AA]=adjacency(DD);
%ii=find((ismember(DD.types,{'u','r','f','c','p'}) & any(AA(:,clist)>0,2)') | ...
%         ismember(DD.types,{'f'}));
ii=false(1,size(AA,1)); 
ii(clist)=true;
ii = ii | any(AA(:,clist)>0,2)';
ii(plist)=false;
ii=find(ii);

parents=getparents(DD);
ni=length(ii);
F=cell(1,ni);
Fvars=cell(1,ni);
for i=1:ni
  var=ii(i);
  cpd=DD.cpds{var};
  if any(strcmp(cpd.type,'v')) || (any(strcmp(cpd.type,'f')) && ~getfunc)
    F{i}=DD.values{var};
    Fvars{i}=fliplr(parents{var});
  elseif isfield(cpd,'cpt')
    if ~isfield(cpd,'order') || any(strfind(cpd.order,'l'))
      Fvars{i}=fliplr(parents{var});
    else
      Fvars{i}=parents{var};
    end
    if getfunc || ~ismember(var,clist)
      F{i}=cpd.cpt;
      if ~isfield(cpd,'order') || any(strfind(cpd.order,'c'))
        Fvars{i}=[var Fvars{i}];
      else
        Fvars{i}=[Fvars{i} var];
      end
    else
      F{i}=DD.values{var}'*cpd.cpt;
    end
  else
    error('should not happen')
  end
  if issparse(F{i}) && forcefull==1
    F{i}=full(F{i});
  end
end
n=DD.sizes;
if getfunc
  F=[{[]} F];
  Fvars=[{fliplr(clist)} Fvars];
end
[Ef,V,cost,Iexpand]=sumproduct(F,Fvars,n,[],fliplr(plist),options);
if isnumeric(Ef)
  Ef=Ef(:);
end

% pass factor names to orderdisp for display purposes
if abs(orderdisplay)>=1 && ~isempty(V)
  if orderdisplay<0, tex=true;  open='^{'; close='}';
  else               tex=false; open='';   close='';
  end
  names=cell(1,length(ii)+1);
  names{1}='V';
  for i=1:length(ii)
    names{i+1}=['P' open num2str(ii(i)) close];
  end
  s=orderdisp(V,names,tex);
  disp(s)
end





