% getenterexit Determines when noise variables enter and exit an EV processing sequence
% USAGE
%   [eenter,eexit]=getenterexit(parents,e);
% INPUTS
%   parents  : m-element cell array with lists of integers for each state variable: 
%                positive integers represent X values 
%                negative integers represent e values
%   e        : de-element cell array of rv structures (discrete or 
%                w/ discrete approximations)
% OUTPUTS
%   eenter   : vector indicating the first variable to use e(i) 
%   eexit    : vector indicating the last variable to use e(i) 
function [eenter,eexit] = getenterexit(parents,e)
de = length(e);
m=length(parents);
Ae = false(m,de);
for i=1:m
  Ae(i,-parents{i}(parents{i}<0))=true;   % Ae(i,j)=1 if e variable j conditions target variable i
end
eenter = zeros(de,1);   % iteration where e(i) enters combined parents
eexit  = zeros(de,1);   % iteration where e(i) exits combined parents
for j=1:de
  ii=find(Ae(:,j),1,'first');
  if ~isempty(ii), eenter(j) = ii;  end
  ii  = find(Ae(:,j),1,'last'); 
  if ~isempty(ii), eexit(j) = ii;  end
end 
