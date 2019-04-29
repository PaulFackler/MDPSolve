% Kumarmom computes the mean and CV of the Kumaraswamy distribution
% USAGE
%   [mu,CV]=Kumar(a,b);
% 
% This is the inverse of Kumarmom2param
% see: Kumarmom2param
function [mu,CV,m2]=Kumarmom(a,b)

mu  = b.*beta(1+1./a,b);
mu2 = b.*beta(1+2./a,b);
CV=sqrt(mu2./mu^2-1);
m2=[mu mu2];

