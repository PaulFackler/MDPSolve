% multioutput Allows anonymous functions to have multiple outputs
% Define an anonymous function as
%   f = @(x) multioutput{f1,f2,...,fd},x);
% and call using
%   [y1,y2,...,yd] = f(x);
% This will make yi = fi(x);
% d must be at least as great as nargout

function [varargout] = multioutput(f,varargin)
varargout = cell(1,nargout);
for i=1:nargout
  varargout{i} = f{i}(varargin{:});
end
