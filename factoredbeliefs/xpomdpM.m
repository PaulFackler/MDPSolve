% xpomdpM Sets up a marginalized extended POMDP 
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
%   EZO   : conditional expected value of the observed variables in Z
%
% The ordering of Xb are first the observable variable, then the actions, 
% then the belief states. 

% This function replaces the unobservable states with belief states  
% measuring the belief weight that a particular value of the unobservable
% state is currently true.

function [Pb,Rb,Sb,Xb,Ixb,EF]=xpomdpM(inc,P,R,X,indox,indux,Z,indoz,induz,passive,F)
if nargin<10, passive=0; end
inda=1:size(X,2); inda(:,[indox indux])=[];   % columns of X with actions
[KJ,nox,nu]=getindex(X,[inda indox],indux);
indy=1:size(Z,2); indy(:,[indoz induz])=[];   % columns of Z with signals
[KI,noz,nu2]=getindex(Z,[indoz indy],induz);

if nu~=nu2
  error('the possible unobservable variables must be the same in X and Z')
end

[Iy,Y]=getI(Z,indy); 
ny=size(Y,1);
[Ios,OS]=getI(Z,indoz); %#ok<ASGLU>
nos=size(OS,1);
if nos*ny~=noz
  %error('not implemented for models in which the observable state and signal are not a regular grid')
end

[IZO,ZO]=getI(Z,[indoz indy]);
KY=getindex(ZO,1:length(indoz),length(indoz)+1:size(ZO,2));

d=length(induz); 
q=zeros(1,d);
for i=1:d, q(i)=length(unique(Z(:,induz(i)))); end
Bf=factoredBgrid(q,inc);
B=B2S(Bf,q,inc);
nb=size(B,1);
nxb=nb*nox;
nzb=nb*noz;
nsb=nb*nos;
% compute the updated belief weights
Bplus=zeros(nos,nxb,nu);
Pb=sparse(nsb,nxb);
EF=0;
radd=repmat(0:nb:nb*(nos-1),nu,nxb);
cadd=repmat(1:nxb,nu*nos,1);
for k=1:ny
  for i=1:nos
    for j=1:nu
      bp=kroncol(reshape(P(KI(KY(i,k),j),KJ),nox,nu),B);
      bp=sum(reshape(bp,nxb,nu),2);
      Bplus(i,:,j)=bp;
      clear bp;
    end
  end
  Bplus=reshape(Bplus,nos*nxb,nu);
  omega=sum(Bplus,2);
  Bplus=Bplus./(omega*ones(1,nu));
  ind=omega==0;
  if any(ind)
    Bplus(ind,:)=0;
  end
  if passive
    Bplus=zeros(nos,nxb,nu);
    for i=1:nos
      Bplus(i,:,:)=repmat(B,nox,1);
    end
    Bplus=reshape(Bplus,nos*nxb,nu);
  end
  Pb=Pb+reshape(mxv(factoredBbas(S2B(Bplus,q,inc),q,inc),omega),nsb,nxb);
  Bplus=reshape(Bplus,nos,nxb,nu);
  Bplus(:)=0;
  if nargout>=10
    ef=sum(P(Iy==k,:).*F(Iy==k,:),1);
    Ef=Ef+reshape(omega,nos,nxb)'*ef;
  end
end

if nargout>1
  Rb=B*R(KJ)'; 
  Rb=Rb(:);
  [Ixb,Xb]=getI(X,[inda indox]); %#ok<ASGLU>
  Xb=rectgrid(Xb,Bf);
  [Ixb,Sb]=getI(Xb,[length(inda)+1:length(inda)+length(indox) length(indox)+length(inda)+(1:size(Bf,2))]);
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
for i=1:no
  for j=1:nu
    ij=find(Io==i & Iu==j);
    switch length(ij)
      case 1
        K(i,j)=ij;
      case 0
      otherwise
        error('Z must contain unique rows')
    end
  end
end

