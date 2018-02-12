% xpomdpamd Creates an xPOMDP model for a standard AMDP model
% USAGE
%   [Pb,Rb,Sb,Xb,Ixb,Pp,P,X,ox,ux,Z,oz,uz]=xpomdpamd(inc,P,R,X,svars);
% INPUTS
%   standard inputs for AMD model - see documentation for amdp model
% OUTPUTS
%   Pb,Rb,Sb,Xb,Ixb   : standard outputs for AMDP and xPOMDP models
%   Pp                : transition matrix for passive adaptive
%   P,X,ox,ux,Z,oz,uz : inputs for xPOMDP models
%                         these can be passed to xpomdpsim for simulation
% Example:
%   [Pb,Rb,Sb,Xb,Ixb,Pp,P,X,indox,indux,Z,indoz,induz]=xpomdpamd(inc,P,R,X,svars);
%   model=struct('P',Pb,'R',Rb,'d',delta,'Ix',Ixb,'T',T);
%   results=mdpsolve(model);
%   Ixopt=resultsa.Ixopt; 
%   [Xvals,EX]=xpomdpsim(o0,u0,b0,TS,P,X,indox,indux,Z,indoz,induz,Xb,Ixopt); 

function [Pb,Rb,Sb,Xb,Ixb,Pp,P,X,ox,ux,Z,oz,uz]=xpomdpamd(inc,P,R,X,svars)
ne=length(P);
[ns,nx]=size(P{1});
na=nx/ns;
ds=length(svars);
dx=size(X,2);
da=dx-ds;
[Ix,S]=getI(X,svars);
Z=rectgrid((1:ne)',S);
X=rectgrid((1:ne)',X);
if isnumeric(R)
  R=repmat(R,ne,1);
else
  R=cell2mat(R); R=R(:);
end
ox=1+svars;
ux=1;
oz=2:ds+1;
uz=ux;
P=blkdiag(P{:});
[Pb,Rb,Sb,Xb,Ixb]=xpomdp(inc,P,R,X,ox,ux,Z,oz,uz,0);
if nargout>=5
  Pp=xpomdp(inc,P,zeros(size(X,1),1),X,ox,ux,Z,oz,uz,1);
end