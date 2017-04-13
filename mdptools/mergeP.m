% mergeP Merges 2 transition probability matrices
% USAGE
%   P=mergeP(P1,P2,options);
% or
%   P=mergeP(P1,P2,colstoch,nx);
% INPUTS
%   P1       : ns1 x nx1 transition probability matrix
%   P2       : ns2 x nx2 transition probability matrix
%   options  : a structure variable with possible fields colstoch and nx
%   colstoch : 1 if P1 and P2 are in column stochastic form [default: 1]
%   nx       : combined number of state/action combinations (only needed if
%                P1 or P2 is a function)
% OUTPUT
%   P  : ns1*ns2 x nx1*nx2 column stochastic matrix (regardless of colstoch)
%
% Alternatively either P1 or P2 can be a cell array with either nx2 or nx1 elements.
% If the probability of S2 depends on X1, P2 should be a cell array in which the
%   ith element is an ns2 x nx2 matrix with probability conditioned on X1(i,:)
% If the probability of S1 depends on X2, P1 should be a cell array in which the
%   ith element is an ns1 x nx1 matrix with probability conditioned on X2(i,:)
%
% An additional option is to pass either P1 or P2 as a function handle. This function
%   should accept a single integer argument and return the conditional probability
%   given the ith value of the other variable. Essentially it should return P1{i} 
%   or P2{i} if these are defined as cell arrays. The main advantage to this
%   syntax is that all of the conditional matrices do not need to be stored.
%   If this syntax is used one must specify nx (the total number of state/action 
%   combinations).
%
% This procedure will create transition matrices that are compatible with
% S and X created using
%   S=rectgrid(S1,S2);
% and 
%   X=rectgrid(X1,X2);
% To create transition matrices composed of more than two subcomponents use
% mergeP recrursively; e.g.
%   P=mergeP(P1,P2); P=mergeP(P,P3);

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function P=mergeP(P1,P2,varargin)
% get options: colstoch and nx
if nargin<3
  colstoch=1; 
  nx=[];
elseif isstruct(varargin{1})
  options=varargin{1};
  if isfield(options,'colstoch'), colstoch=options.colstoch;
  else                            colstoch=1;
  end
  if isfield(options,'nx'), nx=options.nx;
  else                      nx=[];
  end
else
  if ~isempty(varargin{1}), colstoch=varargin{1}; else colstoch=1; end
  if length(varargin)>1 && ~isempty(varargin{2}), nx=varargin{2}; 
  else                                            nx=[];
  end
end

if iscell(P1)
  if ~colstoch, P2=P2.'; end
  cols2=size(P2,2);
  if numel(P1)~=cols2
    error('Inputs are not compatible')
  end
  if isempty(nx)
    nx=0;
    for i=1:cols2
      nx=nx+size(P1{i},1+colstoch);
    end
  end
  if issparse(P1{1}) || issparse(P2)
    nn=0;
    for i=1:cols2
      nn=nn+sum(P2(:,i)~=0)*nnz(P1{i});
    end
    P=sparse([],[],[],size(P1{1},2-colstoch)*size(P2,1),nx,nn);
  else
    P=zeros(size(P1{1},2-colstoch)*size(P2,1),nx);
  end
  for i=1:cols2
    k=i:cols2:nx;
    if colstoch
      P(:,k)=kron(P1{i},P2(:,i));
    else
      P(:,k)=kron(P1{i}.',P2(:,i));
    end
  end
elseif iscell(P2)
  if ~colstoch, P1=P1.'; end
  cols1=size(P1,2);
  if numel(P2)~=cols1
    error('Inputs are not compatible')
  end
  if isempty(nx)
    nx=0;
    for i=1:cols1
      nx=nx+size(P2{i},1+colstoch);
    end
  end
  if issparse(P2{1}) || issparse(P1)
    nn=0;
    for i=1:cols1
      nn=nn+sum(P1(:,i)~=0)*nnz(P2{i});
    end
    P=sparse([],[],[],size(P2{1},2-colstoch)*size(P1,1),nx,nn);
  else
    P=zeros(size(P2{1},2-colstoch)*size(P1,1),nx);
  end
  k=0;
  for i=1:cols1
    k=k(end)+(1:size(P2{i},2));
    if colstoch
      P(:,k)=kron(P1(:,i),P2{i});
    else
      P(:,k)=kron(P1(:,i),P2{i}.');
    end
  end
elseif isa(P1,'function_handle')
  if ~colstoch, P2=P2.'; end
  cols2=size(P2,2);
  if isempty(nx)
    error('Must pass nx if P1 or P2 is defined by a function')
  end
  P1i=P1(1);
  if issparse(P1i) || issparse(P2)   % need to approximate the amount of memory needed
    P=sparse([],[],[],size(P1i,2-colstoch)*size(P2,1),nx,nnz(P2)*nnz(P1i));
  else
    P=zeros(size(P1i,2-colstoch)*size(P2,1),nx);
  end
  clear P1i
  for i=1:cols2
    k=i:cols2:nx;
    if colstoch
      P(:,k)=kron(P1(i),P2(:,i));
    else
      P(:,k)=kron(P1(i).',P2(:,i));
    end
  end
elseif isa(P2,'function_handle')
  if ~colstoch, P1=P1.'; end
  cols1=size(P1,2);
  if isempty(nx)
    error('Must pass nx if P1 or P2 is defined by a function')
  end
  P2i=P2(1);
  if issparse(P2i) || issparse(P1)
    P=sparse([],[],[],size(P2i,2-colstoch)*size(P1,1),nx,nnz(P1)*nnz(P2i));
  else
    P=zeros(size(P2i,2-colstoch)*size(P1,1),nx);
  end
  clear P2i
  k=0;
  for i=1:cols1
    k=k(end)+(1:size(P2(i),2));
    if colstoch
      P(:,k)=kron(P1(:,i),P2(i));
    else
      P(:,k)=kron(P1(:,i),P2(i).');
    end
  end
else
  if colstoch
    P=kron(P1,P2);
  else
    P=kron(P1.',P2.');
  end
end

