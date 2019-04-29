% xpomdp Sets up an extended POMDP 
% USAGE
%   [Pb,Rb,Sb,Xb,Ixb,EZO]=xpomdp(inc,P,R,X,indox,indux,Z,indoz,induz);
% INPUTS
%   inc   : number of beliefs intervals for the belief variables
%   P     : n x m probability matrix for future state and signal (Z) 
%             conditioned on the current state and action (X)
%   R     : m-vector of expected utility conditioned on the current state and action (X)
%   X     : m x q matrix of state/action variable values
%   indox : index of observable state variables in X (elements on {1,...,q})
%   indux : index of unobservable state variables in X (elements on {1,...,q})
%   Z     : n x r matrix of state/signal variable values
%   indoz : index of observable state variables in Z (elements on {1,...,r})
%   induz : index of unobservable state variables in Z (elements on {1,...,r})
% OUTPUTS
%   Pb    : belief augmented transition matrix
%   Rb    : belief augmented reward vector
%   Sb    : matrix of state values for augmented model
%   Xb    : matrix of state/action values for augmented model
%   Ixb   : state variable index vector: Sb(Ixb,:)=Xb(:,columns with states)
%
% The ordering of Xb are first the action variables, then the observable 
% state variables, then the belief states. 
%
% This function replaces the unobservable states with belief states  
% measuring the belief weight that a particular value of the unobservable
% state is currently true.

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

function [Pb,Rb,Sb,Xb,Ixb]=xpomdp(inc,P,R,X,indox,indux,Z,indoz,induz,passive)
if nargin<10, passive=0; end
inda=1:size(X,2); inda(:,[indox indux])=[];   % columns of X with actions
[KJ,nox,nu]=getindex(X,[inda indox],indux);
indy=1:size(Z,2); indy(:,[indoz induz])=[];   % columns of Z with signals
[KI,noz,nu2]=getindex(Z,[indoz indy],induz);

if nu~=nu2
  error('the possible unobservable variables must be the same in X and Z')
end

[Iy,Y]=getI(Z,indy); %#ok<ASGLU>
ny=size(Y,1);
[Ios,OS]=getI(Z,indoz); %#ok<ASGLU>
nos=size(OS,1);
if nos*ny~=noz
  %error('not implemented for models in which the observable state and signal are not a regular grid')
end

[IZO,ZO]=getI(Z,[indoz indy]); %#ok<ASGLU>
KY=getindex(ZO,1:length(indoz),length(indoz)+1:size(ZO,2));

B=simplexgrid(nu,inc,1);
nb=size(B,1);
P=full(P);
% compute the updated belief weights
Pb=sparse(nb,nos*nb*nox);
for k=1:ny
  Bplus=tprodm(B,[2 -5;nb nu],P(KI(KY(:,k),:),KJ),[1 4 3 -5;nos nu nox nu],[nos*nb*nox,nu]);
  omega=sum(Bplus,2);
  if passive  % B does not change
    Bplus=repmat(reshape(B,[1 nb 1 nu]),[nos,1,nox,1]);
    Bplus=reshape(Bplus,nos*nb*nox,nu);
  else
    Bplus=bsxfun(@rdivide,Bplus,omega);
    Bplus(omega==0,:)=0;
  end
  temp=simplexbas(Bplus,nu,inc,1,0);
  temp=mxv(temp,omega,0);
  clear omega
  Pb=Pb+temp;
  clear temp
  if ny>1, fprintf('.'); end
end
if ny>1, fprintf('\n'); end
clear Bplus omega
Pb=reshape(Pb,nos*nb,nox*nb);
Rb=B*R(KJ)'; 
Rb=Rb(:);
if nargout>1
  [Ixb,Xb]=getI(X,[inda indox]); %#ok<ASGLU>
  Xb=rectgrid(Xb,B);
  [Ixb,Sb]=getI(Xb,[length(inda)+1:length(inda)+length(indox) length(indox)+length(inda)+(1:nu)]);
end



% Returns an (no x nu) index matrix K where K(i,j) is the
% row of X associated with the ith value of the observable
% variables (in columns indo of X) and the jth value of 
% the unobservable variables (in columns indu of X).
function [K,no,nu]=getindex(X,indo,indu)
[Io,O]=getI(X,indo);
no=size(O,1);
clear O
[Iu,U]=getI(X,indu);
nu=size(U,1);
clear U
K=zeros(no,nu);
K(Io+no*(Iu-1))=1:no*nu;
if any(K==0)
  error('X must contain unique rows')
end
K=reshape(K,no,nu);