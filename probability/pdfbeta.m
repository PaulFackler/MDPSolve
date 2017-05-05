function p=pdfbeta(x,a,b)
p=x.^(a-1).*(1-x).^(b-1)./beta(a,b);