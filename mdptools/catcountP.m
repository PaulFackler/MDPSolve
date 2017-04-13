% catcountP Creates probability matrices for site count models
% Computes exact probabilities (use catcountPapprox for normal approximation)
% USAGE
%   [P,S,X]=catcountP(N,n,m,p,X,outputtype,pattern);  
% INPUTS
%   N  : number of items
%   n  : number of categories
%   m  : number of category/treatment combinations
%   p  : nxm probability matrix (non-negative values with columns summing to 1)
%           Alternatively p can be a function handle that accepts an
%             m-vector of non-negative integers that sum to N (rows of X)
%             If p has two input arguments then pass Z as well
%   X  : qxm matrix of state/action combinations (each row must sum to N)
%          If omitted all possible combinations are used.
%  outputtype : 0) dense matrix
%               1) sparse matrix
%               2) q-element cell array of R-element column vectors
%               3) q-element cell array of R-element sparse column vectors
%  pattern    : if p is a function and returns a (effectively) block diagonal 
%                 matrix, pass in an n x m matrix with the block diagonal
%                 structure to use the block diagonal solver (this can
%                 result in dramatic speedups).
%  Z          : q-row matrix. When p is a 2-input function X(i,:) and Z(i,:) 
%                 are passed to p
% OUTPUTS
%   P  : Rxq state transition matrix
%   S  : Rxn matrix of grid values
%   X  : qxm matrix of state/action combinations
% where
%   R=(N+n-1)!/N!(n-1)!
% The rows of P are determined by a lexicographic ordering of the points in
%  {0,...,N}^n which sum to N.
% The columns in X correspond to the columns of p.

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2014, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function [P,S,X]=catcountP(N,n,m,p,X,outputtype,pattern,Z)
  if nargin<8, Z=[]; end          
  if nargin<7, pattern=[]; end    
  if nargin<6, outputtype=[]; end  
  if nargin<5, X=[]; end
  if isempty(outputtype), outputtype=0; end % default is dense matrix
  if isnumeric(p)
    pvariable=false;
    p=full(p);
    [mp,np]=size(p);
    if mp~=n,  error('p must have n rows'); end
    if np~=m,  error('p must have m columns'); end
    if ~ischar(pattern)
      if isempty(pattern) 
        blocksr=blkdiagcheck(p);
      else
        blocksr=blkdiagcheck(pattern);
      end
      if length(blocksr)>1
        P=catcountBlkDiag(N,n,m,p,X,outputtype,[]); 
        if nargout>=2, S=simplexgrid(n,N,N,1); end
        if nargout>=3, X=simplexgrid(m,N,N,1,'uint32'); end
        return;
      end
    end
  elseif ~isempty(pattern)
    [mp,np]=size(pattern);
    if mp~=n,  error('pattern must have n rows'); end
    if np~=m,  error('pattern must have m columns'); end
    P=catcountBlkDiag(N,n,m,p,X,outputtype,pattern,Z);
    if nargout>=2, S=simplexgrid(n,N,N,1); end
    if nargout>=3, X=simplexgrid(m,N,N,1,'uint32')'; end
    return
  else    
    pvariable=true;
    if isempty(Z), useZ=false; else useZ=true; end
  end
  
  S=simplexgrid(n,N,N,0,'uint32')';
  % X is not input, get all possible values of X
  if isempty(X)
    if m~=n
      X=simplexgrid(m,N,N,0,'uint32');
    else
      X=S';
    end
    X=[X N-sum(X,2)];
    q=size(X,1);
  else
    [q,nm]=size(X);
    if nm~=m,  error('X must have m columns'); end
  end
  
  % if MEX file doesn't exist use the multisum procedure
  if exist('getpzc','file')~=3 
    P=multisum(p,X,outputtype,Z);
    if nargout>1, S=S'; S=[S N-sum(S,2)]; end
    return;
  end
  
  % indices and numbers of the points that sum to i or less (i=0:N)
  xind=cell(N+1,1);
  nxind=uint32(zeros(N+1,1));
  for i=0:N
    xind{i+1} =uint32(find(sum(S,1)<=i))-1;
    nxind(i+1)=uint32(length(xind{i+1}));
  end
  
  tab=gettab(N,max(n,m));         % used to determine indices
  factor=cumsum(log([1 1:N]))';   % used in computing probabilities
  wl0=warning('off','MATLAB:log:logOfZero');
  if ~pvariable
    logp=log(p);             
  end
  switch outputtype         % initialize the output
    case 0
      P=zeros(size(S,2),q);
    case 1
      P=sparse(size(S,2),q);        % initialize the output
    case {2,3}
      P=cell(1,q);                  % initialize the output
    otherwise
      error('invalid outputtype')
  end
  % loop over the q rows of X
  for j=1:q
    Xj=X(j,:)';
    if pvariable
      if useZ
        logp=log(full(p(double(Xj),Z(j,:)')));
      else
        logp=log(full(p(double(Xj))));
      end
    end
    Pj=getpzc(S,xind,nxind,tab,factor,logp,uint32(Xj));  % MEX file for speed
    switch outputtype
      case 0
        P(:,j)=Pj;
      case 1
        P=add2sparse(P,Pj,j,[],1);
        %P(:,j)=Pj;
      case 2 
        P{j}=Pj;
      case 3
        P{j}=sparse(Pj);
    end
  end
  warning(wl0)
  if nargout>1, S=S';  S=double([S N-sum(S,2)]); end


% table of factorial values
% Used in computing the indices
function tab=gettab(N,m)
  tab=zeros(N+2,m);
  tab(:,1)=(0:N+1)';
  for j=2:m
    tab(:,j)=cumsum(tab(:,j-1));
  end
  if tab(end,end)>intmax('uint32')
    error('problem size is too large')
  end
  tab=uint32(tab);
  
% Specialized version of catcountP to process block diagonal individual
% site transition matrices
function P=catcountBlkDiag(N,n,m,p,X,outputtype,pattern,Z)
S=simplexgrid(n,N,N);
if isempty(X)
  X=simplexgrid(m,N,N);
end
en=eye(n);
un=ones(n,1);
if isnumeric(p)
  [blocksr,blocksc]=blkdiagcheck(p);
  p=full(p);
  nblocks=length(blocksr);
  em=eye(m);
  um=ones(m,1);
  for i=1:nblocks
    br=blocksr{i};
    bc=blocksc{i};
    ni=length(br);
    mi=length(bc);
    pi=[p(br,bc) zeros(ni,1);zeros(1,mi) 1];
    Pi=catcountP(N,ni+1,mi+1,pi,[],0,'skip');
    ss=un; ss(br)=0; ss=[en(:,br)  ss]; %#ok<*AGROW>
    ind=simplexindex(S*ss,ni+1,N,N);
    Pi=Pi(ind,:);
    if outputtype==1
      Pi=sparse(Pi);
    end
    ss=um; ss(bc)=0; ss=[em(:,bc)  ss];
    ind=simplexindex(X*ss,mi+1,N,N);
    if i>1
      Pi=Pi(:,ind);
      P=P.*Pi;
    else
      P=Pi(:,ind);
    end
    clear Pi
  end
else
  q=size(X,1);
  switch outputtype
    case 0
      P=zeros(size(S,1),q);
    case 1
      P=sparse(size(S,1),q);
    case {2,3}
      P=cell(1,q);
  end
  [blocksr,blocksc]=blkdiagcheck(pattern);
  nblocks=length(blocksr);
  ind=cell(1,nblocks);
  Si=cell(1,nblocks);
  xind=cell(1,nblocks);
  nxind=cell(1,nblocks);
  for i=1:nblocks
    br=blocksr{i};
    ni=length(br);
    ss=un; ss(br)=0; ss=[en(:,br)  ss]; %#ok<*AGROW>
    ind{i}=uint32(simplexindex(S*ss,ni+1,N,N));
    Si{i}=simplexgrid(ni+1,N,N,0,'uint32')';
    % indices and numbers of the points that sum to i or less (i=0:N)
    xind{i}=cell(N+1,1);
    nxind{i}=uint32(zeros(N+1,1));
    for j=0:N
      xind{i}{j+1} =uint32(find(sum(Si{i},1)<=j))-1;
      nxind{i}(j+1)=uint32(length(xind{i}{j+1}));
    end
  end
  tab=gettab(N,max(n,m));         % used to determine indices
  factor=cumsum(log([1 1:N]))';   % used in computing probabilities
  wl0=warning('off','MATLAB:log:logOfZero');
  Pj=zeros(size(S,1),nblocks);
  for j=1:q
    Xj=double(X(j,:)');
    if isempty(Z)
      pj=full(p(Xj));
    else
      pj=full(p(Xj,Z(j,:)'));
    end
    for i=1:nblocks
      br=blocksr{i};
      bc=blocksc{i};
      ni=length(br);
      mi=length(bc);
      pi=[pj(br,bc) zeros(ni,1);zeros(1,mi) 1];
      Xi=Xj(br);
      Xi=uint32([Xi;N-sum(Xi)]);
      Pi=getpzc(Si{i},xind{i},nxind{i},tab,factor,log(pi),Xi);  % MEX file for speed
      Pi=Pi(ind{i});
      if i==1
        Pj=Pi;
      else
        Pj=Pj.*Pi;
      end
    end
    if abs(sum(Pj)-1)>1e-14
      error('shouldn''t happen')
    end
    switch outputtype
      case 0
        P(:,j)=Pj;
      case 1
        P=add2sparse(P,Pj,j,[],1);
      case 2 
        P{j}=Pj;
      case 3
        P{j}=sparse(Pj);
    end
  end
  warning(wl0);
end
