% obsupdate constructs a belief augmented transition matrix 
% USAGE
%   Pb=obsupdate(inc,Q,S,indu);
% INPUTS
%   inc  : number of increments used to define the belief states
%   Q    : ny x ns observation probability matrix
%   S    : ns x q matrix of state/action variable values
%   indu : index of columns of S associated with unobservable states
% OUTPUTS
%   Pb  : belief augmented transition matrix

function Pb=obsupdate2(inc,Q,S,indu)
[ny,ns]=size(Q);
if size(S,1)~=ns
  error('# of columns in Q must equal # of rows in S')
end
nv=size(S,2);
nuv=length(indu);
% handle the perfectly observed case - signal doesn't change beliefs
if nuv==0
  Pb=speye(ns);
  return
end
nov=nv-nuv;
indo=1:nv; indo(indu)=[];
 
[Ios,OS]=getI(S,indo);  % indexes the observable states
[Ius,US]=getI(S,indu);  % indexes the unobservable states
nos=size(OS,1);
nus=size(US,1);
clear US
if nus==1
  B=1;
else
  B=simplexgrid(nus,inc,1);
end
nb=size(B,1);
B=rectgrid(OS,B);
clear OS
nsb=size(B,1);
Ix=getI(B,1:nov);
B=B(:,nov+1:end);

% indexes mapping from belief augmented states to original states
KI=invind([Ios Ius]);
KI=KI(Ix,:);
Bplus=zeros(nsb,ny,nus);
omega=zeros(nsb,ny);
if 0
bp=zeros(nus,1);
for i=1:nsb      
  for j=1:ny
    for k=1:nus
      bp(k)=Q(j,KI(i,k))*B(i,k);
    end
    omega(i,j)=sum(bp);
    bp=bp/omega(i,j);
    Bplus(i,j,:)=bp;
  end
end
else
for j=1:ny
  bp=reshape(Q(j,KI),nsb,nus).*B;
  omega(:,j)=sum(bp,2);
  Bplus(:,j,:)=vxm(1./omega(:,j),bp);
  %Bplus(:,j,:)=bp./omega(:,j+zeros(1,nus));
  clear bp
end
end
Bplus(isnan(Bplus))=0;
Pb=sparse(nb,nsb);
for k=1:ny
  Pb=Pb+mxv(simplexbas(reshape(Bplus(:,k,:),nsb,nus),nus,inc,1),omega(:,k));
end
Pb=kroncol(kron(speye(nos),ones(1,nsb/nos)),Pb);

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
