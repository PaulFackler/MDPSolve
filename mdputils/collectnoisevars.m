% collectnoisevars Collects a set of noise variables in e and w 
% A set of variables defined using rvdef are combined into a cell 
% array of noise variables e and a cell array of probability vectors w.
% These can then be passed to g2P.
% USAGE
%   [e,w]=collectnoisevars(v1,v2,...);
% INPUTS
%   v1, v2,..., vd : a set of random variable structures defined using rvdef
% OUTPUTS
%   e   : 1 x d cell array containing vectors of values of the noise variables
%   w   : 1 x d cell array containing associated probability vectors

function [e,w]=collectnoisevars(varargin)
d=length(varargin);
e=cell(1,d);
w=cell(1,d);
for i=1:d
  e{i}=varargin{i}.values;
  w{i}=varargin{i}.cpt;
end
