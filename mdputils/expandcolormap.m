% expandcolormap expands the default colormap to obtain up to 256 color gradations
% USAGE
%   expandcolormap
% The new colormap will interpolate between the old colors to provide a
%   smoother looking plot
function expandcolormap(maxC)
if nargin<1, maxC=256; end
C=colormap;
k=floor(maxC/size(C,1));
if k<2, return; end
C1=C(1:end-1,:);
C2=C(2:end,:);
v=linspace(0,k-1,k)'/k;
C=[1-v v]*[C1(:)';C2(:)'];
C=[reshape(C,k*size(C1,1),3);C2(end,:)];
colormap(C)