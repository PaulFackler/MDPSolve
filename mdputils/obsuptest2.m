
inc=4;     % number of belief subintervals
nu=1;
S=[1;2;3];
ns=size(S,1);  
ny=2;

Q=rand(ny,ns); Q=mxv(Q,1./sum(Q));

%Q=rand(ny,1); Q=Q/sum(Q); Q=Q*ones(1,ns);

%Q=eye(ns,ns);
Q=[0 1 1;1 0 0];
tic
[Pb2]=obsupdate(inc,Q,S,1);
toc

tic
[Pb2b]=obsupdate2(inc,Q,S,1);
toc

max(max(abs(Pb2-Pb2b)))

P=rand(ns,ns); P=mxv(P,1./sum(P));

%P=[1 0 0;1/3 1/3 1/3;1/4 1/2 1/4]';
%P=eye(ns);
tic
[Pb1]=timeupdate(inc,P,S,[],1);
toc
tic
[Pb1b]=timeupdate2(inc,P,S,[],1);
toc
max(max(abs(Pb1-Pb1b)))

R=ones(ns,1);
Qtype=1;
[b,Pb]=pomdp(inc,P,Q,R,struct('Qtype',Qtype));

if Qtype
PP=Pb2*Pb1;
else
PP=Pb1*Pb2;
end
 
 disp(' ')
  max(max(abs(Pb-PP)))
  
  figure(1);spy(PP,'r');hold on; spy(Pb,'b'); hold off
  
  return
  
  PP=Pb2*Pb1;
  figure(2); clf
  [X,Y]=rectgrid((1:size(b,1))',(1:size(b,1))');
  patchplot(X,Y,Pb(:)-PP(:),[-1 1]);
  
 return 

show([zeros(ns,ns) b';b Pb],3)
disp(' ')
show([zeros(ns,ns) b';b Pb2*Pb1],3)
