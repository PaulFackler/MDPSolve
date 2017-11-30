function results = multilogit(y,x,beta0,maxit,tol,print)
% PURPOSE: implements multinomial logistic regression
% Pr(y_i=j) = exp(x_i'beta_j)/sum_l[exp(x_i'beta_l)]
%   where:
%   i    =   1,2,...,nobs
%   j,l  = 0,1,2,...,ncat
%-------------------------------------------------------------------------%
% USAGE: results = logit(y,x,beta0,maxit,tol,print)
% where: y = response variable vector (nobs x 1)
%            the response variable should be coded sequentially from 0 to
%            ncat, i.e., y in {0,1,2,...,ncat}
%        x = matrix of covariates (nobs x nvar)
%            NOTE: to include a constant term in each beta_j,
%            include a column of ones in x
%    beta0 = optional starting values for beta (nvar x ncat+1) (default=0)
%    maxit = optional maximum number of iterations (default=100)
%      tol = optional convergence tolerance (default=1e-6)
%    print = 0/1 1 to print results of each iteration (default=0)
%-------------------------------------------------------------------------%
% RETURNS: a structure
%        results.meth = 'multilogit'
%    results.beta_mat = (nvar x ncat) matrix of beta coefficients:
%                       [beta_1 beta_2 ... beta_ncat] under the 
%                       normalization beta_0 = 0
%    results.beta_vec = (nvar*ncat x 1) vector of beta coefficients:
%                       [beta_1 ; beta_2 ; ... ; beta_ncat] under
%                       normalization beta_0 = 0
%        results.covb = (nvar*ncat x nvar*ncat) covariance matrix
%                       of results.beta_vec
%   results.tstat_mat = matrix of t-statistics conformable to 
%                       results.beta_mat
% results.tstat_vec   = vector of t-statistics conformable to
%                       results.beta_vec
%        results.yfit = (nobs x ncat+1) matrix of fitted 
%                       probabilities: [P_0 P_1 ... P_ncat]
%                       where P_j = [P_1j ; P_2j ; ... ; P_nobsj]
%         results.lik = unrestricted log likelihood
%        results.cnvg = convergence criterion
%        results.iter = number of iterations
%        results.nobs = number of observations
%        results.nvar = number of variables
%        results.ncat = number of categories of dependent variable
%                       (including the reference category j = 0)
%       results.count = vector of counts of each value taken by y, i.e., 
%                       count = [#y=0 #y=1 ... #y=ncat]
%           results.y = y vector
%      results.lratio = LR test statistic against intercept-only model (all 
%                       betas=0), distributed chi-squared with (nvar-1)*ncat
%                       degrees of freedom
%        results.rsqr = McFadden pseudo-R^2
%      
%-------------------------------------------------------------------------%
% A NOTE: Since users might prefer results (coefficients and tstats) in 
%   either a vector or matrix format, and since there is no single natural
%   representation for these in the multinomial logit model, the results
%   structure returns both.  Note that the input arguments require that 
%   (optional) starting values in matrix (nvar x ncat) format.
%
%-------------------------------------------------------------------------%
% SEE ALSO: prt_multilogit, multilogit_lik
%-------------------------------------------------------------------------%
% References: Greene (1997), p.914

% written by:
% Simon D. Woodcock
% CISER / Economics
% 201 Caldwell Hall
% Cornell University
% Ithaca, NY 14850
% sdw9@cornell.edu
%
% modified Nov. 2017, Paul L. Fackler, paul_fackler@ncsu.edu

%---------------------------------------------------------%
%       ERROR CHECKING AND PRELIMINARY CALCULATIONS       %
%---------------------------------------------------------%

if nargin < 2, error('multilogit: wrong # of input arguments'); end;
y = round(y(:)); nobs=size(y,1); [rx, nvar]=size(x);

if (rx~=nobs), error('multilogit: row dimensions of x and y must agree'); end;

% initial calculations
xstd = [1 std(x(:,2:nvar))];
x = x ./ ( ones(nobs,1)*xstd );                             % standardize x
ymin = min(y);
ymax = max(y);
ncat = ymax - ymin;
d0 = ( y*ones(1,ncat+1) ) == ( ones(nobs,1)*(ymin:ymax) );  % put y in dummy format
d = d0(:,2:ncat+1);                                         % normalize beta_0 = 0

% starting values

if nargin < 3
    beta0 = zeros(nvar,ncat+1);
else 
    a = size(beta0,1);
    if a == 0
        beta0 = zeros(nvar,ncat+1);
    else for j = 1:ncat;
            beta0(:,j) = beta0(:,j) .* xstd';
         end;
    end;
end;

beta = beta0(:,2:ncat+1);

% default max iterations and tolerance
if nargin < 4 , maxit = 100; tol = 1e-6; end;
if nargin < 5 , tol = 1e-6; end;
if nargin < 6 , print=0;    end;

if nargin > 7 , error('multilogit: wrong # of arguments'); end;

% check nvar and ncat are consistently defined;
[rbeta, cbeta] = size(beta);
if nvar ~= rbeta
    error('multilogit: rows of beta and columns of x do not agree')
end;
if ncat ~= cbeta
    error(['multilogit: number of columns in beta and categories in y do not agree. ' ...
        'check that y is numbered continuously, i.e., y takes values in {0,1,2,3,4,5}' ...
        ' is ok, y takes values in {0,1,2,3,4,99} is not.'])
end;

%----------------------------------------------------%
% MAXIMUM LIKELIHOOD ESTIMATION OF MULTINOMIAL LOGIT %
%----------------------------------------------------%

% index for block diagonal elements of Hessian
vec=@(x)x(:);
bdind=vec(repmat(reshape(1:nvar*ncat,nvar,ncat),nvar,1)) + ...
      vec(repmat(nvar*ncat*(0:nvar*ncat-1),nvar,1));
d_0 = (y == min(y));

% newton-raphson update
iter=0;
while iter < maxit
    iter = iter + 1;
    [P,lnL] = multilogit_lik(y,x,beta,d,d_0);                   % update P, lnL
    [g, H]  = multilogit_deriv(x,d,P,nobs,nvar,ncat,bdind);     % update g,H
    Hg = H\g(:);
    if abs(g(:)'*Hg/numel(g)) > tol
      beta = beta - reshape(Hg,nvar,ncat);
      if print
        disp(['iteration: ' num2str(iter)]);
        disp(['log-likelihood: ' num2str(lnL)]);
      end
    else
      break; 
    end
end;

%---------------------------------------------------------%
%               GENERATE RESULTS STRUCTURE                %
%---------------------------------------------------------%

results.meth = 'multilogit';
for j = 1:ncat
    results.beta_mat(:,j) = beta(:,j) ./ xstd';        % restore original scale
end;
for j = 1:ncat
    f = (j-1)*nvar + 1;
    l = j*nvar;
    results.beta_vec(f:l,1) = results.beta_mat(:,j);
end;
warning('off','MATLAB:illConditionedMatrix')
results.covb = -inv(H)./kron(ones(ncat),(xstd'*xstd)); % restore original scale
warning('on','MATLAB:illConditionedMatrix')
stdb = sqrt(diag(results.covb));
results.tstat_vec = results.beta_vec./stdb;
for j = 1:ncat                                  
    f = (j-1)*nvar + 1;
    l = j*nvar;
    results.tstat_mat(:,j) = results.tstat_vec(f:l,1);
end;
P_0 = ones(nobs,1) - sum(P,2);
results.yfit = [P_0 P];
results.lik = lnL;
results.cnvg = tol;
results.iter = iter;
results.nobs = nobs;
results.nvar = nvar;
results.ncat = ncat;
results.count = [nobs-sum(sum(d),2) sum(d)];
results.y = y;

% basic specification testing;
p = results.count / nobs;
lnLr = nobs*sum((p.*log(p)),2); % restricted log-likelihood: intercepts only
results.lratio = -2*(lnLr - results.lik);
results.rsqr = 1 - (results.lik / lnLr); % McFadden pseudo-R^2


function [P,lnL] = multilogit_lik(y,x,beta,d,d_0)
% PURPOSE: Computes value of log likelihood function for multinomial logit regression
e_xb = exp(x*beta);
sum_e_xb1 = 1 + sum(e_xb,2);
P=bsxfun(@rdivide,e_xb,sum_e_xb1);
P_0 = 1 - sum(P,2);
lnL = sum(d .* log(P)) + sum(d_0.*log(P_0));


%---------------------------------------------------------%
%     SUPPLEMENTARY FUNCTION FOR COMPUTING DERIVATIVES    %
%---------------------------------------------------------%
function [g,H] = multilogit_deriv(x,d,P,nobs,nvar,ncat,bdind)
% PURPOSE: Computes gradient and Hessian of multinomial logit model
% -----------------------------------------------------------------
% References: Greene (1997), p.914

% compute gradient matrix (nvar x ncat)
g=x'*(d-P);

% compute Hessian, which has (ncat)^2 blocks of size (nvar x nvar)
xp=bsxfun(@times,x,reshape(P,[nobs 1 ncat]));
xp=reshape(xp,[nobs nvar*ncat]); 
H=xp'*xp;
xp=x'*xp;
% subtract off the ncat nvar x nvar blocks from the diagonal blocks
H(bdind)=H(bdind)-xp(:);



% alternative that is slower but possibly more transparent
function [g,H] = multilogit_deriv2(x,d,P,nobs,nvar,ncat,bdind)
% PURPOSE: Computes gradient and Hessian of multinomial logit model
% -----------------------------------------------------------------
% References: Greene (1997), p.914

% compute gradient matrix (nvar x ncat)
g=x'*(d-P);

% compute Hessian, which has (ncat)^2 blocks of size (nvar x nvar)
% this algorithm builds each block individually, m&n are block indices
pp=cell(1,ncat);
for m=1:ncat,
  pp{m}=bsxfun(@times,x,P(:,m)); 
end
H=cell(ncat,ncat);
for m = 1:ncat; 
    H{m,m} = pp{m}'*(pp{m}-x)  ;
    for n = m+1:ncat;
      H{m,n} = pp{m}'*pp{n};
      H{n,m} = H{m,n}'; 
    end;
end;
H=cell2mat(H);