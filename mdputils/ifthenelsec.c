/*
% ifthenelse One line if/then/else function
% USAGE
%   x=ifthenelse(condition,trueres,falseres);
% Inputs
%   condition : logical array
%   trueres   : vector with x(i)=trueres(i) when condition(i)=true
%   falseres  : vector with x(i)=falseres(i) when condition(i)=false
% OUTPUT
%   x         : array of the same size as condition
%
% Note: trueres and falseres can be scalrs or arrays of equal size as
% condition
%
% Performs like the C expression x = c ? t : f
%
% Example:
%   a=randn(3,2,2); b=2; x=ifthenelse(a<b,a,b);
% Returns the same value as x=min(a,b);
*/

#include "mex.h"
#include <math.h>

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[]) {
  double  *Y, *Z, *A, yval, zval;
  mwSize n, i;
  bool *X, yscalar=false, zscalar=false;
  
  /* Error checking on inputs */  
  if (nrhs>3 || nrhs<3) mexErrMsgTxt("Invalid number of input arguments");
  for (i=0; i<nrhs; i++){ 
    if (!mxIsDouble(prhs[i]))
      if (i>0) mexErrMsgTxt("Inputs 2 and 3 must be double");  
    if (mxIsSparse(prhs[i]))
      mexErrMsgTxt("Inputs cannot be sparse");
    if (mxIsComplex(prhs[i]))
      mexErrMsgTxt("Inputs cannot be complex");
  }
  
  n = mxGetNumberOfElements(prhs[0]);
  if (mxGetNumberOfElements(prhs[1]) != n)
    if (mxGetNumberOfElements(prhs[1]) == 1)
      yscalar = true;
    else
      mexErrMsgTxt("X and Y are incompatible");
  if (mxGetNumberOfElements(prhs[2]) != n)
    if (mxGetNumberOfElements(prhs[2]) == 1)
      zscalar = true;
    else
      mexErrMsgTxt("X and Z are incompatible");
   
  X = mxGetData(prhs[0]);
  Y = mxGetPr(prhs[1]);
  Z = mxGetPr(prhs[2]);
  
  plhs[0]=mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
  A=mxGetPr(plhs[0]);

  if (zscalar){
    zval = *Z;
    if (yscalar){
      yval = *Y;
      for (i=0; i<n; i++) A[i] =  X[i] ? yval : zval; 
    }
    else {
      for (i=0; i<n; i++) A[i] =  X[i] ? Y[i] : zval; 
    }
  }
  else {
    if (yscalar){
      yval = *Y;
      for (i=0; i<n; i++) A[i] =  X[i] ? yval : Z[i]; 
    }
    else {
      for (i=0; i<n; i++) A[i] =  X[i] ? Y[i] : Z[i]; 
    }
  }
}


