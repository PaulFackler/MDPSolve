#include "mex.h"
#include <math.h>
/*
% samplenoreplacec utility for samplenoreplace 
% USAGE
%   I=samplenoreplace(e,m);
% INPUTS
%   e  : p x n matrix of uniform [0,1] random values
%   m  : # of elements in set
% OUTPUT
%   I  : p x n mattrix of values on {1,...,m}
% Each column contains p unique values 
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double  *e, *I;
  mwSize m, n, p, i, j, i1, k, ii, ik, m1, *ind;

  /* Error checking on inputs */  
  if (nrhs!=2) mexErrMsgTxt("Not enough input arguments");
  for (i=0; i<nrhs; i++) {
    if (!mxIsDouble(prhs[i]) && !mxIsSparse(prhs[i]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[i]))
      mexErrMsgTxt("e  must be real.");
  }
  
  p=mxGetM(prhs[0]);
  n=mxGetN(prhs[0]);
  e=mxGetPr(prhs[0]);
  m=(mwSize)*mxGetPr(prhs[1]);
  if (p>m)
      mexErrMsgTxt("rows(e) must be <= m");
      
  plhs[0]=mxCreateDoubleMatrix(p,n,mxREAL);
  I = mxGetPr(plhs[0]);
  ind=mxMalloc(m*sizeof(mwSize));
  for (i=0; i<m; i++) ind[i]=i+1;
  m1=m-1;
  for (j=0; j<n; j++){
    for (i=0; i<p; i++) {
      k=floor(*e++ * (m-i));
      if (k>=m || k<0) {
        printf("%1i\n",k);
        mexErrMsgTxt("invalid k");
      }
      ik=ind[k];
      ind[k]=ind[m1-i];
      ind[m1-i]=ik;
      *I++ = (double) ik;
    }
  }
  mxFree(ind);
}


