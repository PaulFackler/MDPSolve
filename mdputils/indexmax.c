#include "mex.h"
#include <math.h>
/*
% indexmax Determines the maximum relative to an index
% USAGE
%   [vmax,mind] = indexmax(v,ind,n);
% INPUTS
%   v   : m-vector
%   ind : m-vector of index values between 1 and n
%   n   : scalar positive integer
% OUTPUT
%   vmax : n-vector of values
%   mind : n-vector of indices 
%
% vmax(i)=max(v(ind==i));
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double  *v, *vplus, *sind, vi;
  mwSize q, n, i, si;
  int ii;
  unsigned int *x;

  /* Error checking on inputs */  
  if (nrhs!=3) mexErrMsgTxt("Not enough input arguments");
  for (ii=0; ii<nrhs; ii++) {
    if (!mxIsDouble(prhs[ii]) && !mxIsSparse(prhs[ii]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[ii]))
      mexErrMsgTxt("X must be real.");
  }
  
  q=mxGetNumberOfElements(prhs[0]);
  if (mxGetNumberOfElements(prhs[1])!=q)
      mexErrMsgTxt("Inputs must have the same number of elements");

  v   =mxGetPr(prhs[0]);
  v--;
  sind=mxGetPr(prhs[1]);
  sind--;
  n   =*mxGetPr(prhs[2]);
      
  plhs[0]=mxCreateDoubleMatrix(n,1,mxREAL);
  vplus=mxGetPr(plhs[0]);
  vplus--;
  vi=-mxGetInf();
  for (i=1; i<=n; i++) vplus[i]=vi;
    
  if (nlhs<2)
    for (i=1; i<=q; i++){
      si=sind[i];
      if (v[i]>vplus[si]) {vplus[si]=v[i];}
    }
  else{
    plhs[1]=mxCreateNumericMatrix(n,1,mxUINT32_CLASS,mxREAL);
    x=mxGetData(plhs[1]);
    //plhs[1]=mxCreateDoubleMatrix(n,1,mxREAL);
    //x=mxGetPr(plhs[1]);
    x--;
    for (i=1; i<=q; i++){
      si=sind[i];
      if (v[i]>vplus[si]) {vplus[si]=v[i]; x[si]=i;}
    }
    //x++;
  }
  //mvplus++;
}


