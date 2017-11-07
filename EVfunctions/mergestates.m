% mergestates Merges groups of states
% USAGE
%   [S,P,X,Ix]=mergestates(C1,C2,C3,...,options);
% INPUTS
%   Cj : cell arrays containing some or all of the elements
%           {Sj,Pj,Aj,INDj}
%        Alternatively Cj can be a matrix which is interpreted as
%           Sj, with Pj, Aj and INDj treated as empty
%   options : a structure variable with the following allowed fields:
%               colstoch: 0/1 1 if transition matrices are in colstoch form
% OUTPUTS
%   Merged values of S, P, X and Ix
% where
%   S   : an ns-row matrix of state values
%   P   : an ns x nx transition probability matrix (or nxm if options.colstoch=1)
%   X   : an nx-row matrix of state/action values
%   Ix  : an nx-vector of state indices
%
% This function preserves lexicographic ordering of the states. Specifically
% if the Sj matrices have lexicographically ordered rows, S will also.
%
% See documentation for extended functionality

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2017, Paul L. Fackler (paul_fackler@ncsu.edu)
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
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR X PARTICULAR PURPOSE 
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

function [S,P,X,Ix]=mergestates(varargin)
if isstruct(varargin{end})
  options=varargin{end};
  ngroup=nargin-1;
else
  options=[];
  ngroup=nargin;
end
colstoch=false;
if ~isempty(options)
  if isfield(options,'colstoch'), colstoch=options.colstoch; end
end
[S,P,X,Ix]=unpack(varargin{1},1);
if isa(P,'function_handle') && ngroup>2
  error('Only two groups can be merged if either P is a function')
end
for j=2:ngroup
  [s,p,a,ind]=unpack(varargin{j},j);
  if isa(p,'function_handle') && ngroup>2
    error('Only two groups can be merged if either P is a function')
  end
  [S,P,X,Ix]=merge(S,P,X,Ix,s,p,a,ind,colstoch);
  clear s p a ind
end


function [S,P,X,Ix]=merge(S,P,X,Ix,s,p,a,ind,colstoch)
  N=size(S,1);
  n=size(s,1);
  M=size(Ix,1);
  m=size(ind,1);
  R=kron(speye(N),ones(n,1));
  if isa(P,'function_handle')
    if isa(p,'function_handle')
      error('One P must be a matrix')
    else
      if isempty(a), aa=s; else aa=a; end
      P=condmerge2(aa,p,P,colstoch);
    end
  else
    if isa(p,'function_handle')
      if isempty(X), AA=S; else AA=X; end
      P=condmerge(AA,P,p,colstoch);
    else
      P=kron(P,p);
    end
  end
  S=[R*S repmat(s,N,1)];
  if isempty(X)
    if ~isempty(a), X=repmat(a,N,1); end
  else
    R=kron(speye(M),ones(n,1));
    if ~isempty(a), X=[R*X repmat(a,N,1)];
    else            X=R*X; 
    end
  end
  if isempty(Ix)
    Ix=ind;
  elseif ~isempty(ind)
    Ix=ones(m,1)*((Ix(:)'-1)*n) + ind(:)*ones(1,M);
  end
  Ix=Ix(:);
  
  
% gets the s, p, a and ind elements of the jth group  
% if X is empty Ix(i)=i, i=1,...,n (n is the # of states)
function [S,P,X,Ix]=unpack(cellj,j)
if isnumeric(cellj)
  S=cellj;
  P=[];X=[];Ix=(1:size(S,1))';
else
  switch numel(cellj)  % note that 3 inputs are not allowed
    case 1
      S=cellj{1}; P=[]; X=[]; Ix=(1:size(S,1))';
    case 2
      S=cellj{1}; P=cellj{2}; X=[]; Ix=(1:size(S,1))';
    case 4
      S=cellj{1}; P=cellj{2}; X=cellj{3}; Ix=cellj{4};
    otherwise
      error(['Cell array for group ' num2str(j) ' should contain 1, 2 or 4 elements'])
  end
  Ix=Ix(:);
  if ~isempty(X) && numel(Ix)~=size(X,1)
    error(['X and Ix for group ' num2str(j) ' are incompatible in size'])
  end
end


% condmerge faciliates conditional merging of states
% USAGE
%   P=condmerge(X1,P1,P2,colstoch);
% INPUTS
%   X1 : m1-row matrix state/action combinations for group 1
%   P1 : m1-row matrix of transition probabilities for group 1
%   P2 : function handle to compute group 2 probabilities
%           P2(X1(i,:)') should return a transition matrix for group 2
%   colstoch : 0/1 if 1 all probabilities are assumed to be colstoch
%                   so P1 is an m1-column matrix
% OUTPUT
%   P  : the merged transition probability matrix 
function P=condmerge(X1,P1,P2,colstoch)
  if nargin<3 || nargin>4 
    error('either 3 or 4 inputs must be passed')
  end
  if nargin<4
    colstoch=false;
  end
  if ~isempty(X1) && ~isnumeric(X1) 
    error('X1 must be a matrix')
  end
  if ~isnumeric(P1)
    error('P1 must be a matrix')
  end
  if ~isa(P2,'function_handle') || nargin(P2)~=1
    error('P2 must be defined as a function of a single argument')
  end
  if colstoch
    m1=size(X1,1);
    if m1~=size(P1,2)
      error('If colstoch rows(X1) must equal cols(P1)')
    end
    P=kron(P1(:,1),P2(X1(1,:)'));
    % handle sparse differently from dense P
    if issparse(P)
      for i=2:m1
        P=[P kron(P1(:,i),P2(X1(i,:)'))];  %#ok<AGROW>
      end
    else
      temp=P;
      P=zeros([size(temp) m1]);
      P(:,:,1)=temp;
      clear temp;
      for i=2:m1
        P(:,:,i)=kron(P1(:,i),P2(X1(i,:)'));
      end
      P=reshape(P,size(P,1),size(P,2)*m1);
    end
  else
    m1=size(X1,1);
    if m1~=size(P1,1)
      error('If not colstoch X1 and P1 must have the same number of rows')
    end
    P=kron(P1(1,:),P2(X1(1,:)'));
    % handle sparse differently from dense P
    if issparse(P)
      for i=2:m1
        P=[P;kron(P1(i,:),P2(X1(i,:)'))]; %#ok<AGROW>
      end
    else
      temp=P;
      P=zeros([size(temp,1) m1 size(temp,2)]);
      P(:,1,:)=temp;
      clear temp;
      for i=2:m1
        P(:,i,:)=kron(P1(i,:),P2(X1(i,:)'));
      end
      P=reshape(P,size(P,1)*m1,size(P,3));
    end
  end
  
  
  % condmerge2 faciliates conditional merging of states
% USAGE
%   P=condmerge2(X2,P2,P1,colstoch);
% INPUTS
%   X2 : m1-column matrix state/action combinations for group 2
%   P2 : m1-column matrix of transition probabilities for group 2
%   P1 : function handle to compute group 1 probabilities
%           P1(X1(i,:)') should return a transition matrix for group 1
%   colstoch : 0/1 if 1 all probabilities are assumed to be colstoch
%                   and X2 and P2 are m2-column matrices
% OUTPUT
%   P  : the merged transition probability matrix 
%
% This merges two groups with the first group conditional on the second
% To merge 2 groups with the second conditional on the first use condmerge
function P=condmerge2(X2,P2,P1,colstoch)
  % error('not working correctly')
  if nargin<3 || nargin>4 
    error('either 3 or 4 inputs must be passed')
  end
  if nargin<4
    colstoch=false;
  end
  if ~isnumeric(X2) || ~isnumeric(P2)
    error('The first two inputs must be matrices')
  end
  if ~isa(P1,'function_handle') || nargin(P1)~=1
    error('P1 must be defined as a function of a single argument')
  end
  if colstoch
    m2=size(X2,1);
    if m2~=size(P2,2)
      error('If colstoch rows(X2) must equal cols(P2)')
    end
    Pi=kron(P1(X2(1,:)'),P2(:,1));
    [mi,ni]=size(Pi);
    if issparse(Pi), P=sparse(mi,m2*ni);
    else             P= zeros(mi,m2*ni);
    end
    ind=0:m2:m2*(ni-1);
    P(:,ind+1)=Pi;
    clear Pi
    for i=2:m2
      Pi=kron(P1(X2(i,:)'),P2(:,i));
      P(:,ind+i)=Pi;
      clear Pi
    end
  else
    [m2]=size(X2,1);
    if m2~=size(P2,1)
      error('If not colstoch X2 and P2 must have the same number of rows')
    end
    Pi=kron(P1(X2(1,:)'),P2(1,:));
    [mi,ni]=size(Pi);
    if issparse(Pi),  P=sparse(mi*m2,ni);
    else              P=zeros(mi*m2,ni);
    end
    ind=0:m2:m2*(mi-1);
    P(ind+1,:)=Pi;
    clear Pi
    for i=2:size(X2,1)
      Pi=kron(P1(X2(i,:)'),P2(i,:));
      P(ind+i,:)=Pi; 
      clear Pi
    end
  end