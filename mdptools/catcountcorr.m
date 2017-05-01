p1=[0.2 0.9;0.8 0.1]; p2=[0.8 0.2;0.2 0.8]; w=[0.4;0.6]; p=p1*w(1)+p2*w(2);
N=5;
P1=catcountP(N,2,2,p1);
P2=catcountP(N,2,2,p2);
P=catcountP(N,2,2,p);

Pw=P1*w(1)+P2*w(2);

figure(1); clf
for i=0:N
subplot(1,N+1,i+1)
plot((0:N)',[P(:,i+1) Pw(:,i+1)],'-*')
xlim([0 N])
end
pl=longrunP(P);
plw=longrunP(Pw);
figure(2); clf
plot((0:N)',[pl plw])
