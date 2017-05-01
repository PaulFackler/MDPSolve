
inc=10;     % number of belief subintervals
no=5;
nu=3;
S=rectgrid((1:no)',(1:nu)');
ns=size(S,1);  
ny=5;

P=rand(ns,ns); P=mxv(P,1./sum(P));
P=eye(ns);
P=kron(rand(no,no),eye(nu)); P=mxv(P,1./sum(P));
Q=rand(ny,ns); Q=mxv(Q,1./sum(Q));
%Q=rand(ny,1); Q=Q/sum(Q); Q=Q*ones(1,ns);
%Q=eye(ns,ns);


tic
[Pb1]=timeupdate2(inc,P,S,[1],[2]);
toc
tic
[Pb1b]=timeupdate(inc,P,S,[1],[2]);
toc
max(max(abs(Pb1-Pb1b)))
return

clc
tic
[Pb2b]=obsupdate2(Q,S,[1 2],inc);
toc
tic
[Pb2]=obsupdate(Q,S,[1 2],inc);
toc
max(max(abs(Pb2-Pb2b)))

nnz(Pb1*Pb2), nnz(Pb1)+nnz(Pb2)
return
R=ones(ns,1);
[b,Pb]=pomdp(inc,P,Q);
 
nnz(Pb), nnz(Pb1*Pb2), nnz(Pb1), nnz(Pb2)
return
 disp(' ')
  max(max(abs(Pb-Pb2*Pb1)))
  max(max(abs(Pb-Pb1*Pb2)))
  
  PP=Pb2*Pb1;
  figure(2); clf
  [X,Y]=rectgrid((1:126)',(1:126)');
  patchplot(X,Y,Pb(:)-PP(:),[-1 1]);