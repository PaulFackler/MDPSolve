% obsupdate constructs a belief augmented transition matrix 
% USAGE
%   [Pb]=obsupdate(inc,Q,S,indu);
% INPUTS
%   inc  : number of increments used to define the belief states
%   Q    : ny x ns observation probability matrix
%   S    : ns x q matrix of state/action variable values
%   indu : index of columns of S associated with unobservable states
% OUTPUTS
%   Pb  : belief augmented transition matrix

function [Pb]=obsupdate(inc,Q,S,indu)
[ny,ns]=size(Q);   % # of signal values by # of state values
if size(S,1)~=ns
  error('# of columns in Q must equal # of rows in S')
end
nv=size(S,2);      % # of state variables
nuv=length(indu);  % # of unobservable state variables
% handle the perfectly observed case - signal doesn't change beliefs
if nuv==0,  Pb=speye(ns); return; end
indo=1:nv; indo(indu)=[];
[Ios,OS]=getI(S,indo);  % indexes the observable states
[Ius,US]=getI(S,indu);  % indexes the unobservable states
nos=size(OS,1);
nus=size(US,1);
clear OS US
if nus==1, B=1; else B=simplexgrid(nus,inc,1); end
nb=size(B,1);  % # of belief values
nsb=nb*nos;    % # of augmented state values
% indexes mapping from belief augmented states to original states
KI=invind([Ios Ius]);
Pb=sparse(nb,nsb);
% loop over the values of the signal
for k=1:ny
  bp=kroncol(reshape(Q(k,KI),nos,nus),B);  % updated belief weights
  omega=sum(bp,2);
  ind=omega==0;    % possible that signal k cannot be obtained for some states
  omega(ind)=1;    % reset omega to an arbitrary value to avoid 0/0 problems
  bp=vxm(1./omega,bp);
  %bp=bp./omega(:,ones(1,nus));
  omega(ind)=0;
  Pb=Pb+mxv(simplexbas(bp,nus,inc,1),omega);
  clear bp omega
end
% expand Pb to be block diagonal with 1 block per observable state value
Pb=kroncol(kron(speye(nos),ones(1,nb)),Pb);

% I is an n x 2 matrix with the index of the observable variable in the
% first column and the index of the unobservable in the second.
% There are no and nu values of the observable and unobservable variables
% and K is an no x nu matrix that maps (Oi,Uj) into the index number of the
% combined variables.
function K=invind(I)
n=max(I);
K=zeros(n);
for i=1:n(1)
  for j=1:n(2)
    ij=find(I(:,1)==i & I(:,2)==j);
    switch length(ij)
      case 1
        K(i,j)=ij;
      case 0
      otherwise
        error('I must contain unique rows')
    end
  end
end


