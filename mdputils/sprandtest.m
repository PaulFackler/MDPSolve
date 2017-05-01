N=(10:10:100)';
n=length(N);
U=250;

p=ones(n,U); 
pi=ones(n,1);
p(:,1)=pi./N;
for i=2:U,
  pi=pi+(N-pi)./N;
  p(:,i)=pi./N; 
end; 

figure(1); plot(p,(1:U))