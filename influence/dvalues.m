% dvalues Variable values from an influence diagram
% USAGE
%   X=dvalues(D,list,options);
% INPUTS
%   D       : an influence diagram structure
%   list    : an m-element cell array of variable names or 
%             an m-element vector of variable numbers
%   options : a structure of options
%               matrix      - set to 1 to have output be a matrix
%               passforward - set to 1 to have functions evaluated at
%                                parent values
%               dorder      - lexicographic with the same order as D
%             a string can also be passed
%               'm' is the same as struct('matrix',1)
%               'p' is the same as struct('passforward',1)
%               'd' is the same as struct('dorder',1)
%               'mp' or 'pm' are the same as struct(,'matrix',1,'passforward',1)
% OUTPUTS
%   X    : an m-element cell array or m-column matrix of variable values
%   ind  : index vector of variables involved in the expansion
%             The size of the vectors in X equals prod(D.sizes(ind))
%   DD   : a diagram structure with the variables defined by functions
%             expanded (mainly for internal purposes)
%
% Note that passforward can be set to an array of logicals to allow some 
%   variable functions to be passed forward and others to not be
%
% If any variables are passed forward the ordering of the variables is
%   lexicographic with respect to the ordering of the variables in the
%   diagram (i.e., if the passforward option is used then the dorder option
%   is automatically used regardless of the value of dorder).
%
% To obtain a cell array with the parent values for variable i use:
%   X=dvalues(D,D.parents{i});
% To obtain this as a matrix use:
%   X=dvalues(D,D.parents{i},'m');

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

function [X,ind,D]=dvalues(D,list,options)
dD=length(D.sizes);
matrix=false;
passforward=false(1,dD);
dorder=false;
if nargin>=3 && ~isempty(options)
  if isstruct(options)
    if isfield(options,'matrix'),      matrix=options.matrix;           end
    if isfield(options,'passforward'), passforward=options.passforward; end
    if isfield(options,'dorder'),      dorder=options.dorder;           end
  elseif ischar(options)
    if any(options=='m'), matrix      = true;       end
    if any(options=='p'), passforward = true(1,dD); end
    if any(options=='d'), dorder      = true;       end
  end
end
if isnumeric(list)
  check=ismember(list,1:dD);
else
  [check,list]=ismember(list,D.names);
end
if ~all(check)
  error('list contains variables not in D')
end
if all(size(passforward)==1)
  passforward=passforward(ones(1,dD));
end
options=struct('matrix',false,'passforward',passforward);

% no passforwards needed
if ~any(passforward(list)) & ~any(cellfun(@isempty,D.values(list)))
  if dorder, [list,ind]=sort(list); end
  if matrix
    X=mrectgrid(D.values{list});
    if dorder, X=X(:,ind); end
  else
    X=crectgrid(D.values{list});
    if dorder, X=X{ind}; end
  end
  ind=list;
  return
end

% evaluate any variables defined by functions and determine
% which variables are involved in the list
if ~all(cellfun(@isnumeric,D.parents))
  D.parents=getparents(D);
end
ind=false(1,dD);
sizes=ones(1,dD);
vars=cell(1,length(list));
for j=1:length(list)
  var=list(j);
  cpd=D.cpds{var};
  if strcmp(cpd.type,'f') && passforward(var) 
    % if function values not yet obtained then evaluate now and update cpd
    if ~isfield(cpd,'expvalues')
      [parents,cpd.expparents,D]=dvalues(D,D.parents{var},options);
      cpd.expvalues=double(D.cpds{var}.valfunc(parents{:}));
      D.cpds{var}=cpd;
      %disp(var)
    end
    vars{j}=cpd.expparents;
    ind(cpd.expparents)=true;
    sizes(var)=length(cpd.expvalues);
  else
    vars{j}=var;
    ind(var)=true;
    sizes(var)=length(cpd.values);
  end
end
% expand out the values for the list variables so they have the same size
ind=find(ind);
if matrix
  X=zeros(prod(sizes),length(list));
  for i=1:length(list)
    var=list(i);
    ir=ones(1,dD); 
    ir(vars{i})=D.sizes(vars{i});
    if strcmp(D.cpds{var}.type,'f') && passforward(var)
      Xi=reshape(D.cpds{var}.expvalues,fliplr(ir));
    else
      Xi=reshape(D.values{var},fliplr(ir));
    end
    ir=ones(1,dD); 
    ir(ind)=D.sizes(ind); 
    ir(vars{i})=1;
    Xi=repmat(Xi,fliplr(ir));
    X(:,i)=Xi(:);
    Xi=[]; %#ok<NASGU>
  end
else
  X=cell(1,length(list));
  for i=1:length(list)
    var=list(i);
    ir=ones(1,dD); 
    ir(vars{i})=D.sizes(vars{i});
    if strcmp(D.cpds{var}.type,'f') && passforward(var)
      X{i}=reshape(D.cpds{var}.expvalues,fliplr(ir));
    else
      X{i}=reshape(D.values{var},fliplr(ir));
    end
    ir=ones(1,dD); 
    ir(ind)=D.sizes(ind); 
    ir(vars{i})=1;
    X{i}=repmat(X{i},fliplr(ir));
    X{i}=X{i}(:);
  end
end


function X=mrectgrid(varargin)
d=nargin;
n=1:d;
for i=1:d
  n(i)=numel(varargin{i});
end
nX=prod(n);
X=zeros(nX,d);
c0=[1 cumprod(n)];
c1=[fliplr(cumprod(fliplr(n))) 1];
for i=1:d
  xi=repmat(varargin{i}(:)',[c1(i+1) 1 c0(i)]);
  X(:,i)=xi(:);
  xi=[];
end

function X=crectgrid(varargin)
d=nargin;
n=1:d;
for i=1:d
  n(i)=numel(varargin{i});
end
X=cell(1,d);
c0=[1 cumprod(n)];
c1=[fliplr(cumprod(fliplr(n))) 1];
for i=1:d
  xi=repmat(varargin{i}(:)',[c1(i+1) 1 c0(i)]);
  X{i}=xi(:);
end
