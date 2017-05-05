function f=pdfcauchy(x,mu,sigma)
f=1./(pi*sigma.*(1+((x-mu)./sigma).^2));