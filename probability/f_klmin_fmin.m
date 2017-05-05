%function [F] = f_klmin_fmin(x,intlnsbar,intlnomsbar)
%_fmin...to run on cluster
function F = f_klmin_fmin(x,intlnsbar,intlnomsbar)
% klmin: Kullback-Leibler divergence minimization
% Specify FOCs for the Kullback-Leibler minimization:
% Called from stoch_surv_VFI.m
 
 
%F = [psi(x(1)) - psi(x(1)+x(2)) - intlnsbar;
%     psi(x(2)) - psi(x(1)+x(2)) - intlnomsbar];
 
if min(x)<0  %penalty for non-feasible levels of x
    Fsolve=[inf;inf];
    Jac=NaN*ones(2,2);
else
    cc = psi(x(1)+x(2));
    Fsolve = [psi(x(1)) - cc - intlnsbar;
              psi(x(2)) - cc - intlnomsbar];
end
F=Fsolve'*Fsolve; %for fminsearch must square the output and combine.