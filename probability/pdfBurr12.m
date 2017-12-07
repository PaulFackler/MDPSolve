function p=pdfBurr12(x,a,b,c)
p=a*c/b*(x/b).^(c-1)./(1+(x/b).^c).^(a+1);