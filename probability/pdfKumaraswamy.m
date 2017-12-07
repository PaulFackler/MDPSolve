function f=pdfKumaraswamy(x,a,b)
f=a*b*x.^(a-1).*(1-x.^a).^(b-1);