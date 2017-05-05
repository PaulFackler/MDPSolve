% pomdpsolve Solves finite horizon POMDPs
% Implements the Monahan/Eagle algorithm
% USAGE
%   [V,A]=pomdpsolve(P,Q,R,delta,T,v0,current);
% INPUTS
%   P        : (ns x ns x na) array of state transition probabilities
%                can be (ns x (ns*na)) with na blocks each of size ns x ns stacked horizontally
%   Q        : (ny x ns x na) array of conditional observation probabilities
%                can be (ny x (ns*na)) with na blocks each of size ny x ns stacked horizontally
%                (currently must be ns x na)
%   R        : (na x ny) or (nx x 1) or (nx x ny) array of rewards (see below)
%   delta   : discount factor
%   T       : time horizon
%   v0      : terminal value function (default=0)
%   current : 0/1, 1 if signal is conditioned on current state
% OUTPUTS
%   V : T-element cell array with the feasible alpha vectors
%            each element is a n x k matrix of alpha vectors
%            (k varies over t)
%   A : T-element cell array with each element a k vector
%            of the actions associated with the alpha vectors
% Note: V{t} and A{t} refer to time t so V{1} and A{1} have T periods to go
%  before the terminal date and V{T} and A{T} have 1 period to go.
%
% The value function and optimal actions at belief state lambda 
%    [vt,at]=max(lambda'*V{t}); at=A{t}(at);
%
% Number of variable values
%   ns S (states)
%   na A (actions)
%   ny Y (obervations)
%   nx X (state/action combinations) - nx=ns*na (same # of actions per state)

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

function [V,A]=pomdpsolve(P,Q,R,delta,T,v0,current,Rtype)

[ns,nx]=size(P);  % P is ns x ns*na
na=nx/ns;
if na~=round(na)
  error('P has improper size');
end
ny=numel(Q)/nx;
Q=reshape(Q,[ny,ns,na]);

if nargin<6 || isempty(v0),  v=zeros(ns,1);
else                         v=v0;
end
if nargin<7 || isempty(current), current=0; end
if nargin<8, Rtype=[]; end

% determine Rtype if needed from the size of R
if isempty(Rtype)
    if all(size(R)==[na,ny])
      if all([na,ny]==[ns,na])
        error('Cannot determine if R is ns x na or na x ny - set options.Rtype')
      end
      Rtype=1;
    elseif all(size(R)==[nx,1]) || all(size(R)==[ns,na])
      Rtype=2;
    elseif all(size(R)==[nx,ny])
      Rtype=3;
    else
      error('The size of R is not consistent with other input data')
    end
end

switch Rtype
  case 1  % R is na x ny
    R=repmat(reshape(R,[1,na,ny]),[ns 1 1]);
    R=sum(reshape(R,nx,ny).*reshape(Q,ny,nx)',2);
  case 2  % R is ns x na
    % nothing to do
  case 3  % R is nx x ny
    R=sum(reshape(R,nx,ny).*reshape(Q,ny,nx)',2);
end
R=R(:);

if exist('dompurge','file')==3
  dompurgemex=true;
else
  dompurgemex=false;
end
%dompurgemex=false;
% W is used by getnewv to update of the solution vectors.
% Q=P(Y|S,A) for current and Q=P(Y|S+,A) for future conditioning.
% Note: P is multipied by delta here to avoid future multiplications.
if current
  W=repmat(reshape(delta*P,ns,ns,na,1),[1 1 1 ny]).* ...
    repmat(reshape(permute(Q,[2 3 1]),[1 ns na ny]),[ns 1 1 1]);
else
  W=repmat(reshape(delta*P,ns,ns,na,1),[1 1 1 ny]).* ...
    repmat(reshape(permute(Q,[2 3 1]),[ns 1 na ny]),[1 ns 1 1]);
end
W=reshape(W,ns,ns*na*ny)';

J=[];
C=-ones(1,ns); % dummy coefficient matrix for LP
nonbasis=(1:ns)';

V=cell(T,1);
A=cell(T,1);
% loop over T periods
for t=T:-1:1
  v=getnewv(v);
  [v,a]=prune(v);
  V{t}=v;
  A{t}=a;
end



% get all possible condidate alpha vectors
function vnew=getnewv(v)
  nv=size(v,2);   % number of alpha vectors in next period
  nJJ=nv^ny;       % number of combinations of new vectors
  % expand J if it is not big enough
  if size(J,2)<nJJ
    J=rectgrid(repmat({(1:nv)'},1,ny));
    [jj,ii]=sort(max(J,[],2));
    J=J(ii,:);
  end
  JJ=J(1:nJJ,:);  % indexes all possible combinations
  vtemp=reshape(W*v,nx,ny,nv);
  vnew=zeros(nx,nJJ);
  for i=1:ny
    vnew=vnew+reshape(vtemp(:,i,JJ(:,i)),nx,nJJ);
  end
  vnew=reshape(R*ones(1,nJJ)+vnew,ns,na*nJJ);
end

% prune first prunes alpha vectors that are dominated by other single alpha
% vectors (dompurge). It then solves, for each j,
%   max_w sum(w) s.t.w>=0, sum(w)<=1 and sum_k (v_k-v_j)*w_k <=0
% for all remaining k. If sum(w)<1, v_j is immediately purged.
function [vnew,anew]=prune(v)
  tol=1e-15;
  nv=size(v,2);
  anew=repmat(1:na,1,nv/na);
  % prune the dominated vectors
  if dompurgemex       % if mex file is available it is much faster
    ind=dompurge(v);
    nv=length(ind);
  else                 % if mex file not available
    ind=1:nv;
    for i=nv:-1:1
      B=v(:,ind)-v(:,i)*ones(1,nv);
      B(:,i)=-1;
      if any(all(B>=0)) 
        ind(i)=[]; nv=nv-1;
      end
    end
  end
  v=v(:,ind);     % keep only non-dominated alpha vectors
  anew=anew(ind); % and associated actions
  % prune remaining set using LP
  if nv>2
    ind=1:nv;
    for i=nv:-1:1
      B=v(:,ind)-v(:,i)*ones(1,nv); B(:,i)=1;
      B=[B' zeros(nv,1);C 0]; B(i,end)=1;
      basis=ns+(1:nv)';
      [x,err]=lpx(B,basis,nonbasis);
      if sum(x)<tol || err>0           % prune out element i
        ind(i)=[]; nv=nv-1; 
      end
    end
    v=v(:,ind);
    anew=anew(ind);
  end
  vnew=v;
end

end



%
% Performs simplex steps with a starting basis B for which the
% basis and nonbasic variables are listed in the the vectors
% BASIS and NONBASIC.  
%
function [x,err]=lpx(B,basis,nonbasis)
maxcount=5000;
% Iterate until convergence @
count=0;
err=0; 
tol=-4*eps;
[n1,m1]=size(B);
n=n1-1;
m=m1-1;
[c,j]=min(B(n1,1:m));
while c<tol && (count<=maxcount);
  count=count+1;
  i=find(B(1:n,j)>-tol);
  r=B(i,m1)./B(i,j);
  i=i(r==min(r));
  if size(i,1)>1                       % in degenerate case
    i=i(ceil(rand(1,1))*size(i,1),:);  % pick randomly
  end
  if B(i,j)==0
    disp(' ')
  end
  pivot=1/B(i,j);
  tempj=B(:,j)*(-pivot);
  tempi=B(i,:)*pivot;
  B=B-B(:,j)*tempi;
  B(i,:)=tempi;
  B(:,j)=tempj;
  B(i,j)=pivot;
  tempi=basis(i);
  basis(i)=nonbasis(j);
  nonbasis(j)=tempi;
  [c,j]=min(B(n1,1:m));
end % while
if count>maxcount, err=3; end  % 'Maximum iterations exceeded'
% extract optimal solution 
[basis,rowind]=sort(basis);
B=B(rowind,end);
x=zeros(n,1);
ii=find(basis<=m);
x(basis(ii))=B(ii);
end
