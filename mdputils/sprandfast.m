function x=sprandfast(m,n,p)
N=ceil(m*n*p);
r=ceil(rand(N,1)*m);
c=sort(ceil(rand(N,1)*n));
x=sparse(r,c,rand(N,1),m,n);
