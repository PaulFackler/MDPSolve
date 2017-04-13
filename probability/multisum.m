% multisum Probability distribution of the sum of multinomials
% USAGE
%   [P,S]=multisum(p,X);
% INPUTS
%   p : nxm probability matrix (p>=0 with sum(p)=ones(1,m))
%       Alternatively p can be a handle for a function that takes 
%         a row of X and returns an nxm probability matrix. This allows
%         p to be conditional on X.
%   X : qxm matrix of non-negative integers that sum to N
% OUTPUT
%   P : q-column matrix of probabilities with sum(P)=ones(1,q)
%   S : n-column matrix of values of the multinomial random variable
% The number of rows in P and S is equal to (N+n-1)!/N!/(n-1)!
%
% P(:,i) is the probability distribution of S=sum(z)
%   where z(j) is Multinomial(p(:,j),X(i,j))

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

function [P,S]=multisum(p,X,outputtype)
if nargin<3, outputtype=0; end
N=sum(X,2);
if any(N~=N(1))
  error('All rows of X must sum to the same value')
end
N=N(1);
% handle special case of N=0
if N==0
  if size(X,2)~=size(p,2)
    error('p and X are not comformable')
  end
  P=1; S=zeros(1,size(p,1));
  return
end
wl0=warning('off','MATLAB:log:logOfZero');  % turn log of 0 warning off
if isnumeric(p)
  lnp=log(p);       
  pvariable=false;
else
  lnp=log(p(double(X(1,:))));
  pvariable=true;
end
lnp(lnp==-inf)=-realmax;
[n,m]=size(lnp);
if size(X,2)~=m
  error('p and X are not comformable')
end
% precompute grids and scaling factors
gln=[0 log(cumprod(1:N))]'; % log of factorial function
SS=cell(1,N);
factor=cell(1,N);
for j=1:N
  SS{j}=double(simplexgrid(n,j,j,1)); 
  % factors used in computing multinomial probabilities
  factor{j} = gln(j+1) - sum(gln(SS{j}+1),2);  
end

C=catcountind(n,N);       % precompute indices
pvals=cell(m,N);
  
M=size(X,1);
switch outputtype
  case 0 
    P=zeros(size(SS{end},1),M);
  case 1
    P=sparse(size(SS{end},1),M);
  case {2,3}
    P=cell(1,M);
end
for j=1:M
  Xj=double(X(j,:));
  if pvariable
    lnp=log(p(Xj));
    lnp(lnp==-inf)=-realmax;
    pvals=cell(m,N);
  end
  [Xj,xind]=sort(Xj,2,'descend');
  cxj=Xj(1);  % cumulative sum of Xj
  xi=xind(1);
  if isempty(pvals{xi,cxj})
    pvals{xi,cxj}=exp(SS{cxj}*lnp(:,xi) + factor{cxj});
  end
  Pj=pvals{xi,cxj}; % multinomial probability
  nz=length(Pj);
  for i=2:length(Xj)
    Xji=Xj(i);
    xi=xind(i);
    if Xji==0, break; end
    x=SS{Xji};
    nx=size(x,1);
    if isempty(pvals{xi,Xji})
      pvals{xi,Xji}=exp(x*lnp(:,xi) + factor{Xji});
    end
    Pj=Pj*(pvals{xi,Xji}.');
    zind=C{cxj,Xji};
    nz=zind(nz,nx);
    %Pj=indexsum(Pj,zind,nz);
    Pj=accumarray(zind(:),Pj(:),[nz,1],[],0,false);
    cxj=cxj+Xji;  
  end
  if outputtype>1
    P{j}=Pj;
  else
    P(:,j)=Pj;
  end
end
warning(wl0) % turn log of 0 warning back on
S=SS{N};


% catcountind Precomputes indices for adding category count vectors
% USAGE
%   C=catcountind(d,N);
% INPUTS
%  d : # of categories
%  N : # of sites
% OUTPUTS
%  C : N-1 x N-1 cell array of matrices
%
% If we have the probabilities associated with n1 sites (p1) and the
%   probabiilties associated with n2 sites (p2) then the elements of
%   p1*p2' are distributed according to C(n1,n2) in the probability
%   vector associated with n1+n2 sites.
% Let pi be the probability vector associated with ni sites starting in
%  category i. Then 
%    Ci=C(1:sum(n(1:i-1)),1:ni);
%    Ri=R(i-1)*pj';
%    mi=multiset(d,sum(n(1:i)));
%    Ri=accumarray(Ci,Ri(:),mi,@sum,0);
%  updates R
function C=catcountind(d,N)
C=cell(N-1,N-1);
s=cell(N-1,1);
ns=zeros(N-1,1);
for i=1:N-1
  s{i}=simplexgrid(d,i,i);
  ns(i)=size(s{i},1);
end
tab=simplexindex([],d,N,N);

for i=1:N-1
  if i+i<=N
    z=reshape(bsxfun(@plus,reshape(s{i},[ns(i),1,d]),reshape(s{i},[1,ns(i),d])),ns(i)*ns(i),d); 
    C{i,i}=uint64(reshape(simplexindex(z,d,i+i,i+i,tab),ns(i),ns(i)));
  end
  for j=i+1:N-i
    z=reshape(bsxfun(@plus,reshape(s{i},[ns(i),1,d]),reshape(s{j},[1,ns(j),d])),ns(i)*ns(j),d); 
    C{i,j}=uint64(reshape(simplexindex(z,d,i+j,i+j,tab),ns(i),ns(j)));
    C{j,i}=C{i,j}';
  end
end
 
