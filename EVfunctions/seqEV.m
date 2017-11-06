% seqEV Helper function for evaluating sequential (staged) EV functions
% USAGE
%   f=seqEV(EV);
% INPUT
%   EV : m-element cell array of EV functions for an m-stage transition
% OUTPUT
%   f  : function handle of the form f(V) or f(V,I) that defines 
%          a complete EV function
% For example if P=P2*P1 and EV1(V)=P1'*V and EV2(V)=P2'*V 
% the f(V)=EV1(EV2(V)) and f(V,I)=EV1(EV2(V),I)
function f=seqEV(EV)
m=length(EV);
f=@(varargin) seqEVeval(EV,m,varargin{:});

function y=seqEVeval(EV,m,varargin)
y=EV{end}(varargin{1});
for j=m-1:-1:2
  y=EV{j}(y);
end
% check if the evaluation is indexed or full
if length(varargin)>1
  y=EV{1}(y,varargin{2});
else
  y=EV{1}(y);
end