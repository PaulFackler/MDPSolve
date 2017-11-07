% ztrans Performs a z-transformation normalizing the mean and variance
% USAGE
%   z=ztrans(y);
function z=ztrans(y)
z=bsxfun(@minus,y,mean(y,1)); 
z=bsxfun(@rdivide,z,sqrt(mean(z.*z,1))); 