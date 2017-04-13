% amdp Converts an MDP to an Adaptive (belief state) MDP
% USAGE
%   [b,Pb,Rb,Sb,Xb,Ixb]=amdp(p,P,R,S,X,Ix);
% INPUTS
%   p     : number of belief states for each state (scalar)
%   P     : q-element cell array of state transition matrices
%   R     : q-element cell array of reward matrices
%   S     : ns-row matrix of state variable values
%   X     : nx-row matrix of state/action variable values
%   Ix    : index-vector mapping state/actions to states
% OUTPUTS
%   b     : belief nodal values [Nxq]
%   Pb    : belief state transition matrix
%   Rb    : belief state reward matrix
%   Sb    : matrix of augmented state variable values
%   Xb    : matrix of augmented state/action variable values
%   Ixb   : augmented index-vector mapping state/actions to states
%
% Either or both P and R are cell arrays with q elements.
% These define q alternative models over which a belief space is defined.
% This procedure expands the state space to include a set of belief grid 
% values as additional states. The state space dimension becomes ns*N, where 
%   N=(p+q-1)!/p!/(q-1)!
%   Pb     : N*ns x N*nx matrix replacing original ns x nx matrices
%   Rb     : N*nx--vector or N*ns x na matrix replacing original matrices
%   S      : N*ns x (q+K) matrix replacing original nsxK matrix
%              If S is omitted it is assumed to equal (1:ns)'
%
% The expanded state variables are ordered with the belief states last.
% For example if the original model has two state variables and there are
% three alternative models the states are ordered S1, S2, B1, B2, B3
% where Bi is the degree of belief in model i and Si is the original state i
% (note that the B values sum to 1). 

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

function [b,Pb,Rb,Sb,Xb,Ixb]=amdp(p,P,R,S,X,Ix)

if ~exist('R','var'),     R    =[]; end
if ~exist('S','var'),     S    =[]; end
if ~exist('X','var'),     X    =[]; end
if ~exist('Ix','var'),    Ix   =[]; end

% Determine number of models and ensure that R and P are compatible
% If both R and P are scalar the model is not adaptive
if isnumeric(P)
  P={P};
end
if iscell(P)
  qp=numel(P);           % number of alternative models for P
  if qp==1, probtype=0;
  else      probtype=1;
  end
else
  error('Incorrect P type')
end

if isnumeric(R)
  R={R};
end
if iscell(R)
  qr=numel(R);           % number of alternative models for R
  if qr>1, probtype=probtype+2; end
else
  error('Incorrect R type')
end

switch probtype
  case 0              % P and R are both singleton
    Pb=P{1};
    Rb=R{1};
  case 1              % R is singleton
    q=qp;
    R=repmat(R,q,1); 
  case 2              % P is singleton
    q=qr;
    P=repmat(P,q,1);
  case 3              % P and R are both model dependent
    if qr~=qp
      error('R and P are incompatible')
    end
    q=qr;
end

% problem is not adaptive
% both R and P are matrices or contain only 1 element
if probtype==0
  if nargout>=4
    if nargin>=4
      Sb=S;
    else
      Sb=(1:size(R,1))';
    end
  end
  warning('AMDP:ProblemNotAdaptive','Problem is not adaptive')
  return 
end

[ns,nx]=size(P{1});
if nx<ns
  warning('AMDP:TransposedOnly','AMDP implemented for transposed P matrices only!!')
end
b=simplexgrid(q,p,1,true);  % discretized belief space grid
if nargout==1, return; end

N=size(b,1);
warning('off','SIMPLEXBAS:Extrapolation')
% There are three subfunctions to compute Pb
% The first is fully vectorized but this may induce out-of-memory errors
% The second loops over beliefs, the third over actions
% These are tried in order and, if one fails due to a memory error
% the next is tried.
try
  %Pb=zeros(20000,20000);  % debug statement: introduce memory error to force catch
  Pb=getPb(q,p,ns,nx,N,P,b); 
catch AMDPError
  try
    if ~strcmp(AMDPError.identifier,'MATLAB:nomem')
      rethrow(AMDPError)
    end
    disp('AMDP: Standard approach to compute P failed; switching to low memory approach - may take some time')
    %Pb=zeros(20000,20000);    % debug statement: introduce memory error to force catch
    Pb=getPbAlt1(q,p,ns,nx,N,P,b); 
  catch AMDPError
    if ~strcmp(AMDPError.identifier,'MATLAB:nomem')
      rethrow(AMDPError)
    end
    disp('AMDP: Low memory approach to compute P failed; switching to very low memory approach - may really take some time')
    Pb=getPbAlt2(q,p,ns,nx,N,P,b);
  end
end
warning('on','SIMPLEXBAS:Extrapolation')
if nargout==2, return; end

if isempty(R{1})
  Rb=[];
else
  Rb=kron(R{1}(:),b(:,1));
  for i=2:q
    Rb=Rb+kron(R{i}(:),b(:,i));
  end
  Rb=reshape(Rb,numel(Rb)/size(R{1},2),size(R{1},2));
end
if nargout==3, return; end

if isempty(S),  S=(1:ns)'; end
Sb=[kron(speye(ns),ones(N,1))*S repmat(b,ns,1)];

if isempty(X),  Xb=[];
else 
  Xb=[kron(speye(nx),ones(N,1))*X kron(ones(nx,1),speye(N))*b];
end

if isempty(Ix),  Ixb=[];
else 
  Ixb=ones(N,1)*((Ix(:)'-1)*N) + ((1:N)')*ones(1,nx);
  Ixb=Ixb(:);
end

% basic (vectorized) prodedure to obtain Pb
function Pb=getPb(q,p,ns,nx,N,P,b)
if q==2  % save some time and memory
  ww   =      kron(P{1},b(:,1)');
  wws  = ww + kron(P{2},b(:,2)');
  ind=find(wws);  % eliminate zero probability transitions
  wws=wws(ind);
  ww=ww(ind)./wws;
else
  ww=cell(1,q-1);
  ww{1}=kron(P{1},b(:,1)');
  wws=ww{1};
  for i=2:q-1
    ww{i}=kron(P{i},b(:,i)');
    wws=wws+ww{i};
  end
  wws=wws+kron(P{q},b(:,q)');
  ind=find(wws);              % eliminate zero probability transitions
  wws=wws(ind);
  for i=1:q-1
    ww{i}=ww{i}(ind)./wws;
  end
  ww=[ww{:}];
end
% need to make ww full because simplexbas doesn't support sparse
ww=full(ww);
[Pb,ii]=simplexbas(ww,q,p,1);        % get interpolation values 
clear ww
wws=full(wws);
Pb=mxv(Pb,wws); 
clear wws
jj=ones(q,1)*ind';
ii=double(ii)+N*rem(jj-1,ns);        % adjust ii and jj to avoid reshaping Pb from NxnNm to NnxNm
jj=ceil(jj/ns); 
Pb=sparse(ii,jj,Pb,ns*N,nx*N);
clear ii jj



% low memory procedure to obtain Pb - loops over belief values
function Pb=getPbAlt1(q,p,ns,nx,N,P,b)
Pb=sparse(ns*N,N*nx);  
startcol=1;
for j=1:N
  if q==2  % save some time and memory
    ww  =  P{1}*b(j,1);
    wws  = P{2}*b(j,2);
    wws=wws+ww;
    ind=find(wws);  % eliminate zero probability transitions
    wws=wws(ind);
    ww=ww(ind);
    ww=ww./wws;               % normalize transition weights
  else
    wws=P{1}*b(j,1);
    if issparse(wws)
      ww=sparse([],[],[],ns*nx,q-1,nnz(wws)*(q-1));
    else
      ww=zeros(ns*nx,q-1);
    end
    ww(:,1)=wws(:);
    for i=2:q-1
      wwi=P{i}*b(j,i);
      ww(:,i)=wwi(:);
      wws=wws+wwi;
      clear wwi
    end
    wws=wws+P{q}*b(j,q);
    ind=find(wws);                % eliminate zero probability transitions
    wws=wws(ind);
    ww=ww(ind,:);
    ww=vxm(1./wws,ww);            % normalize transition weights
  end
  % need to make ww full because simplexbas doessn't support sparse
  ww=full(ww);
  [Pbj,ii]=simplexbas(ww,q,p,1);       % get interpolation values
  clear ww
  wws=full(wws);
  Pbj=mxv(Pbj,wws); % same as Pbj=Pbj(:).*wws(jj(:))';
  clear wws
  jj=ind(ones(q,1)*(1:length(ind)));  % jj refers to the column of the nonzeros, get column # of original
  clear ind
  ii=double(ii)+N*rem(jj-1,ns);  % adjust ii and jj to avoid reshaping Pb from NxnNm to NnxNm
  jj=ceil(jj/ns);    
  Pbj=sparse(ii,jj,Pbj,N*ns,nx); 
  clear ii jj
  Pb=add2sparse(Pb,Pbj,startcol,0.6,true);
  clear Pbj 
  startcol=startcol+nx;
end
Pb=Pb(:,reshape(1:nx*N,nx,N)');


% very low memory prodedure to obtain Pb
% the assumption here is that N<<nx
% loops over state/action combinations
function Pb=getPbAlt2(q,p,ns,nx,N,P,b)
Pb=sparse(ns*N,nx*N);
startcol=1;
for j=1:nx
  if q==2  % save some time and memory
    ww   = reshape(P{1}(:,j)*b(:,1)',ns*N,1);
    wws  = reshape(P{2}(:,j)*b(:,2)',ns*N,1);
    wws=wws+ww;
  else
    wws=P{1}(:,j)*b(:,1)';
    if issparse(wws)
      ww=sparse([],[],[],ns*N,q-1,nnz(wws)*(q-1));
    else
      ww=zeros(ns*N,q-1);
    end
    ww(:,1)=wws(:);
    for i=2:q-1
      wwi=P{i}(:,j)*b(:,i)';
      ww(:,i)=wwi(:);
      wws=wws+wwi;
      clear wwi
    end
    wws=wws+P{q}(:,j)*b(:,q)';
  end
  ind=find(wws);  % eliminate zero probability transitions
  wws=wws(ind);
  ww=ww(ind,:);
  ww=full(vxm(1./wws,ww));            % normalize transition weights
  [Pbj,ii]=simplexbas(ww,q,p,1);  % get interpolation values
  clear ww
  wws=full(wws);
  Pbj=mxv(Pbj,wws);  % same as Pbj=Pbj(:).*wws(jj(:))';
  clear wws
  jj=ind(ones(q,1)*(1:length(ind)));    % jj refers to the column of the nonzeros, get column # of original
  clear ind
  ii=double(ii)+N*rem(jj-1,ns);  % adjust ii and jj to avoid reshaping Pb from NxnNm to NnxNm
  jj=ceil(jj/ns);
  Pbj=sparse(ii,jj,Pbj,N*ns,N);
  clear ii jj
  Pb=add2sparse(Pb,Pbj,startcol,0.6,true);
  clear Pbj
  startcol=startcol+N;
end