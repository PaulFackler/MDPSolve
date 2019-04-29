% Burr12mom computes the mean, CV and skewness of the Burr-12 distribution
% USAGE
%   [mu,CV,gamma]=Burr12mom(a,b,c);
% 
% This is the inverse of Burr12mom2param
% see: Burr12mom2param
function [mu,CV,gamma,m3]=Burr12mom(a,b,c)

mu  =    a.*b.*beta(1+1./c, b-1./c);
mu2 = a.^2.*b.*beta(1+2./c, b-2./c);
mu3 = a.^3.*b.*beta(1+3./c, b-3./c);
CV=sqrt(mu2./mu.^2-1);
gamma = (mu3-mu.^3.*(3*CV.^2+1))./(mu.*CV).^3;
m3=[mu mu2 mu3];

