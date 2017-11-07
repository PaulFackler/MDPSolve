% evalmodel
% INPUTS
%  xelist : m-element cell array with lists of integers for each state variable: 
%             positive integers represent X values 
%             negative integers represent e values
% OUTPUTS
%   type     : 0 indicates all states variables have unique conditioning variables
%              1 indicates that all state variables have unique conditioning shocks
%              2 indicates a general model
%   Ax       : logical matrix with A(i,j)=true if X(j) is associated with S(i)
%   Ae       : logical matrix with A(i,j)=true if e(j) is associated with S(i)
%   shocksum : shock j can be summed out after state shocksum(j) is processed
%                shocksum=0 unless type=2
%   reordere : shocks should be reordered so the first to be summed out are highest in order
function [type,Ax,Ae,shocksum,reordere]=evalmodel(xelist,options)
if isnumeric(xelist), xelist={xelist}; end
m=length(xelist);
order=[];
if exist('options','var') && ~isempty(options)
  if isfield(options,'order'),     order = options.order;     end
end
if ~isempty(order)
  if length(order)~=m
    error(['options.order must be an ' num2str(m) '-vector']);
  end
  xelist=xelist(order);
end
maxx=0;
maxe=0;
for i=1:m
  if max(xelist{i})>maxx, maxx=max(xelist{i}); end
  if min(xelist{i})<maxe, maxe=min(xelist{i}); end
end
maxe=-maxe;
Ax = false(m,maxx);
Ae = false(m,maxe);
for i=1:m
  Ax(i,xelist{i}(xelist{i}>0))=true;
  Ae(i,-xelist{i}(xelist{i}<0))=true;
end
if maxx>0, [blocksrx,blockscx]=blkdiagcheck(Ax); end
if maxe>0, [blocksre,blocksce]=blkdiagcheck(Ae); end
type=2;                    % general type
if maxe==0 || length(blocksre)==m 
  type=1;                  % unique shocks
  if maxx==0 || length(blocksrx)==m 
    type=0;                % unique conditioning variables
  end
end

shocksum=zeros(1,maxe);
for j=1:maxe
  if sum(Ae(:,j))>1, shocksum(j) = find(Ae(:,j),1,'last'); end
end

[~,reordere]=sort(shocksum,'descend');





