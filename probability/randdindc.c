#include "mex.h"
#include <math.h>

// linear search - simple though not efficient for large p
double findx(double u, double *p){
double x=0;
  while (u>=0){  
    x++;
    u -= *p++;  
  }
  return(x);
}

// finds x such that c[x] <= u < c[x+1] 
void  getx(double *x, double *u, double *p, double *index, size_t n, size_t q){
double *pi;  
size_t i;
  p -= n; // to accommodate base-1 index
  for (i=0; i<q; ++i) {
    pi = p + (n * (size_t)index[i]); // pi is a pointer to the top of column index[i]
    x[i] = findx(u[i], pi);
  }
}

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double  *u, *p, *index, *x;
  size_t  q, n;
  int ii;

  /* Error checking on inputs */  
  if (nrhs!=3) mexErrMsgTxt("Not enough input arguments");
  for (ii=0; ii<nrhs; ii++) {
    if (!mxIsDouble(prhs[ii]) && !mxIsSparse(prhs[ii]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[ii]))
      mexErrMsgTxt("X must be real.");
  }
  
  q=mxGetNumberOfElements(prhs[0]);
  n=mxGetM(prhs[1]);
  
  u = mxGetPr(prhs[0]);
  p = mxGetPr(prhs[1]);
  index = mxGetPr(prhs[2]);
       
  plhs[0]=mxCreateDoubleMatrix(q,1,mxREAL);
  x=mxGetPr(plhs[0]);
  
  getx(x,u,p,index,n,q); 
}

