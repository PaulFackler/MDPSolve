function x=randgamma(n,lambda,theta)
x=gammaincinv(rand(n,1),lambda)*theta;