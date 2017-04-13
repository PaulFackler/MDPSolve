% xpomdpsim Simulates states and actions for the xPOMDP model
% USAGE
%   [Xvals,EX]=xpomdpsim(o0,u0,b0,T,P,X,indox,indux,Z,indoz,induz,Xb,Ixopt);
% INPUTS
%   o0    : k row matrix on initial observable state values
%   u0    : k row matrix on initial unobservable state values
%   b0    : k row matrix on initial belief values
%   T     : time horizon
%   P     : ns x nx probability matrix
%   X     : nx row matrix of state/action values
%   indox : index of columns of X containing the observable states
%   indux : index of columns of X containing the unobservable states
%   Z     : nz row matrix of state/signal values
%   indoz : index of columns of Z containing the observable states
%   induz : index of columns of Z containing the unobservable states
%   Xb    : matrix of observable state/action/belief values
%   Ixopt : index of observable state/belief values
% OUTPUTS
%   Xvals : 5 element cell array composed of k x T+1 matrices of index
%             values
%   EX    : 5 element cell array composed of T+1 row matrices of expected
%             values
% The 5 elements are:
%     1) observable states
%     2) unobservable states
%     3) actions
%     4) signals
%     5) beliefs
%
% P,X,indox,indux,Z,indoz,induz are the same variables passed to xpomdp
% Xb is returned from xpomdp
% Ixopt is obtained from the results returned from mdpsolve

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

function [Xvals,EX]=xpomdpsim(o0,u0,b0,T,P,X,indox,indux,Z,indoz,induz,Xb,Ixopt)
rep=size(o0,1);
if rep==0, rep=size(u0,1); end
if rep==0, rep=size(b0,1); end
nz=size(Z,1);

% column indices for actions and signals
inda=1:size(X,2); inda(:,[indox indux])=[];   % columns of X with actions
indy=1:size(Z,2); indy(:,[indoz induz])=[];   % columns of Z with signals

% index values for all original variables
[Iox,O]=getI(X,indox);
[Iux,U]=getI(X,indux);
[Ia,A] =getI(X,inda);

Ioz=getI(Z,indoz);
Iuz=getI(Z,induz);
[Iy,Y]=getI(Z,indy);

% redefine X and Z to be 3 column matrices: X=[O A U] and Z=[O Y U]
[XX,xind]=sortrows([Iox Ia Iux]);
[ZZ,zind]=sortrows([Ioz Iy Iuz]);
P=P(zind,xind);
%Pcum=full(cumsum(P,1));  % cumulative probability

% get indices for beliefs
indbx=length([indox inda])+1:size(Xb,2);
B=simplexgrid(length(indbx),find(Xb(:,end)==0,1,'first')-1,1);
%Ib=repmat((1:size(B,1))',size(Xb,1)/size(B,1),1);
%[Ib,B]=getI(Xb,indbx);
nb=size(B,1);
nu=length(indbx);

treeXX = kdtree(XX);
treeZZ = kdtree(ZZ);
treeB  = kdtree(B);

Aopt=getI(Xb,1:length(inda));
Aopt=Aopt(Ixopt);

% vector of ones for computing random numbers
u=ones(nz,1);
uu=zeros(rep,1);

% get initial values
oval=match(o0,O);
uval=match(u0,U);

yval=NaN;
bval=B(match(b0,B),:);

Xvals=repmat({zeros(rep,T+1)},1,5);
Xvals{5}=zeros(rep,nu,T+1);
for t=1:T+1
  Xvals{1}(:,t)=oval; 
  Xvals{2}(:,t)=uval;
  Xvals{4}(:,t)=yval;
  Xvals{5}(:,:,t)=bval;
  aval=getaction(oval,bval);
  Xvals{3}(:,t)=aval;
  [onew,uval,yval]=updatestates(oval,uval,aval); 
  bval=updatebeliefs(oval,bval,aval,onew,yval); 
  oval=onew;
end

% compute expectations if requested
if nargout>1
  EX=cell(1,5);
  EX{1}=mean(reshape(O(Xvals{1},:),[rep,T+1,length(indox)]),1);
  EX{1}=reshape(EX{1},[T+1,length(indox)]);
  EX{2}=mean(reshape(U(Xvals{2},:),[rep,T+1,length(indux)]),1);
  EX{2}=reshape(EX{2},[T+1,length(indux)]);
  EX{3}=mean(reshape(A(Xvals{3},:),[rep,T+1,length(inda)]),1);
  EX{3}=reshape(EX{3},[T+1,length(inda)]);
  Xvals{4}(:,1)=1;
  EX{4}=mean(reshape(Y(Xvals{4},:),[rep,T+1,length(indy)]),1); 
  EX{4}=reshape(EX{4},[T+1,length(indy)]);
  EX{4}(1,:)=NaN+zeros(1,length(indy));
  EX{5}=mean(Xvals{5},1);
  EX{5}=reshape(EX{5},[length(indbx),T+1])';
end

% use nearest neighbor beliefs to determine the action
function aval=getaction(oval,bval)
  ii=kdtree_closestpoint(treeB, bval);
  aval=Aopt(ii+nb*(oval-1));
end

% generate random values of the states
function [onew,uval,yval]=updatestates(oval,uval,aval)
  ind=kdtree_closestpoint(treeXX, [oval aval uval]);
  r = rand(1,rep); 
  %indnew = 1+sum(u*r>Pcum(:,ind),1);
  indnew=randdiscc(P,r,ind);
  uval=ZZ(indnew,:);
  onew=uval(:,1);
  yval=uval(:,2);
  uval=uval(:,3);
end
  
% update the belief weights
function bnew=updatebeliefs(oval,bval,aval,onew,yval)
  bnew=zeros(rep,nu);
  indz=kdtree_closestpoint(treeZZ, [onew yval uu])-1;
  indx=kdtree_closestpoint(treeXX, [oval aval uu])-2;
  for i=1:nu
    bi=0;
    for j=1:nu
      bi=bi+P(indz+i+nz*(indx+j)).*bval(:,j);
    end
    bnew(:,i)=bi;
  end
  bnew=vxm(1./sum(bnew,2),bnew);
end
   
end
