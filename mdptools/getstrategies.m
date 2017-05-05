% getstrategies Constructs a matrix of possible strategies
% USAGE 
%   s=getstrategies(X,svars);
% INPUTS
%   X : nx x dx matrix of state/action combinations
%   svars : list of columns of X associated with state variables
% OUTPUT
%   s : ns-column strategies matrix
%
% There are ns combinations of values of the state variables
% Each row of s provides an index vector associated with a strategy
% Thus X(s(i,:),:) is an ns x dx matrix displaying the state/actions
% combinations associated with strategy i.

function s=getstrategies(X,svars)
[Ix,S]=getI(X,svars);
ns=size(S,1);
s=cell(1,ns);
ns=1;
for i=1:ns
  s{i}=find(Ix==i);
  ns=ns*length(s{i});
end
try
  s=rectgrid(s{:});
catch
  disp('Could not create strategies matrix')
  disp(['There are ' num2str(ns) ' possible stratgies'])
  disp('Returning cell array of possible actions instead')
end
  