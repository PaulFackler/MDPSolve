% conditional Constructs the conditional distribution for a subset of variables
% USAGE
%   P=conditional(D,clist,plist,options);
% INPUTS
%   D       : an influence diagram structure
%   clist   : a cell array of variable names of output variables 
%                This must be a subset of D.names.
%                Alternatively a list of numbers corresponding to the
%                variables.
%   plist   : a cell array of variable names of input (conditioning) variables
%                This must be a subset of D.names
%                Alternatively a list of numbers corresponding to the
%                variables.
%   options : structure variable with control options (described below)
% OUTPUTS
%   P        : conditional probability matrix
%   V        : (m-1)x7 cell array containing information on the order
%                factors are processed and the variables involved
%                Use with orddisp function.
%   cost     : total processing cost (number of multiply operations)
%   Iexpand  : If requested Iexpand is an index vector of length
%                np=size(dvalues(D,plist,'m'),1) that can be used to expand
%                the columns of P if some of the variables in plist are
%                not ancestors of the clist variables.
%                Thus 
%                  P=P(:,Iexpand);
%                has the proper size. If not requested Iexpand is applied
%                to the P matrix before returning it.
%   DD       : the diagram after preprocessing
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

function [P,V,cost,Iexpand,DD]=conditional(D,clist,plist,options)
% default option values
if nargin<4, options=[]; end
passforward     = 1;   % pass functions forward to children
cleanup         = 0;   % used to handle extrapolation
orderdisplay    = 0;   % displays sum-product order information
orderonly       = 0;   % 1 to return order info (f and Iexpand set to [])
forcefull       = 0;   % 1 to force sumproduct to use full factors
spthreshold     = 0.5; % number on [0,1] converts to sparse if sparsity
                       %   ratio is less than threshold
if nargin>=3 && ~isempty(options)
  if isfield(options,'passforward'),  passforward=options.passforward;   end
  if isfield(options,'cleanup'),      cleanup=options.cleanup;           end
  if isfield(options,'orderdisplay'), orderdisplay=options.orderdisplay; end
  if isfield(options,'orderonly'),    orderonly=options.orderonly;       end
  if isfield(options,'forcefull'),    forcefull=options.forcefull;       end
  if isfield(options,'spthreshold'),  spthreshold=options.spthreshold;   end
end

% future states are default children  
if nargin<2 || isempty(clist)
  clist=find(ismember(D.types,{'f'}));
end
% current states and actions are default parents
if nargin<3 || isempty(plist)
  plist=find(ismember(D.types,{'s','d','a'}));
end
  
if ~isnumeric(clist), [junk,clist]=ismember(clist,D.names); end %#ok<*ASGLU>
if ~isnumeric(plist), [junk,plist]=ismember(plist,D.names); end

if length(passforward)==1
  passforward=passforward(1,ones(1,length(D.sizes)));
end
passforward(clist)=false;

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
    if ~isfield(cpd,'cpt') && ~passforward(var)
      DD=f2cpt(DD,var,cleanup,passforward);
    end
  end
end

parents=getparents(DD);
[A,AA]=adjacency(DD);
%ii=find((ismember(DD.types,{'u','r','f','c','p'}) & any(AA(:,clist)>0,2)') | ...
%         ismember(DD.types,{'f'}));
       
ii=false(1,size(AA,1)); 
ii(clist)=true;
ii = ii | any(AA(:,clist)>0,2)';
ii(plist)=false;
ii=find(ii);
F=DD.cpds(ii);
n=DD.sizes;
Fvars=cell(1,length(ii));
for i=1:length(ii)
  var=ii(i);
  if isfield(cpd,'cpt')
    if ~isfield(cpd,'order') || any(strfind(cpd.order,'l'))
      Fvars{i}=[var fliplr(parents{var})];
    else
      Fvars{i}=[var parents{var}];
    end
  else
    Fvars{i}=[var fliplr(parents{var})];
  end
  if ~isfield(DD.cpds{var},'cpt') || isempty(DD.cpds{var}.cpt)
    DD=f2cpt(DD,var,cleanup,passforward);
    F{i}=DD.cpds{var}.cpt;
  else    
    F{i}=DD.cpds{var}.cpt;
  end
  if issparse(F{i}) && forcefull==1
    F{i}=full(F{i});
  end
end

%%
[P,V,cost,Iexpand]=sumproduct(F,Fvars,n,fliplr(clist),fliplr(plist),options);
% make sparse if this will help
if nnz(P)/numel(P)<spthreshold
   P=sparse(P);
end
if nargout<4 && ~isempty(Iexpand), P=P(:,Iexpand); end

% pass factor names to orderdisp for display purposes
if abs(orderdisplay)>=1
  if orderdisplay<0, tex=true;  open='^{'; close='}';
  else               tex=false; open='';   close='';
  end
  names=cell(1,length(ii));
  for i=1:length(ii)
    names{i}=['P' open num2str(ii(i)) close];
  end
  s=orderdisp(V,names,tex);
  disp(s)
end