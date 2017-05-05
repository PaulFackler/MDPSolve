function f=pdfdexpon(x,mu,sigma)
f=exp(-abs(x-mu)./sigma)./(2*sigma);