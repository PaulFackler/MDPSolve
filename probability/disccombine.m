% disccombine Combines a set of discrete probability distributions
% USAGE
%  [E,W]=disccombine(e,w,varargin);
% INPUTS
%   e : m-element cell array containing vectors with length n(i)
%   w : m-element cell array containing probability vectors with length n(i)
% or
%   e1,w2,e2,w2,... : each pair of inputs is a value/probability vector pair
% OUTPUTS
%   E : prod(n) x m matrix containing all tuples of values
%   W : prod(n) x 1 probability vector associated w/ E
% 
% This procedure assumes that the m noise terms are mutually independent

function [E,W]=disccombine(e,w,varargin)
if nargin>2
  e={e,varargin{1:2:end-1}};
  w={w,varargin{2:2:end}};
end
E=rectgrid(e);
W=prod(rectgrid(w),2);