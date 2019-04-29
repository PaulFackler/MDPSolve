function p=pdfBurr12(x,a,b,c)
p=b*c/a*(x/a).^(c-1)./(1+(x/a).^c).^(b+1);