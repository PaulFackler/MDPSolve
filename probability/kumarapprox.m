function kapprox = kumarapprox
kmom = @(a,b,p) b.*beta(1+p./a,b);

a=linspace(0.1,2,101)';
b=linspace(0.05,0.85,101)';
[A,B]=rectgrid(a,b);

mu=kmom(A,B,1);
CV=sqrt(kmom(A,B,2)./mu.^2 - 1);
X=crossorthmom([mu CV],1);
theta = X\[A B];

kapprox = @(mu,CV) crossorthmom([mu CV],1)*theta;