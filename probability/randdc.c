#include "mex.h"
#include <math.h>

// finds x such that c[x] <= u < c[x+1] 
void  getx(double *x, double *u, double *c, double *index, double m, size_t p){
double *xend;  
  xend = x + p;
  while (x < xend) {
    *x = index[(size_t) ceil(m * *u)];
    if (*u++ > c[(size_t) *x]) *x += 1;
    ++x;
  }
}

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double  *u, *c, *index, m, *x, *xend;
  size_t  p;
  int ii;

  /* Error checking on inputs */  
  if (nrhs!=3) mexErrMsgTxt("Not enough input arguments");
  for (ii=0; ii<nrhs; ii++) {
    if (!mxIsDouble(prhs[ii]) && !mxIsSparse(prhs[ii]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[ii]))
      mexErrMsgTxt("X must be real.");
  }
  
  p=mxGetNumberOfElements(prhs[0]);
  m=(double)mxGetNumberOfElements(prhs[2]);
  
  u = mxGetPr(prhs[0]);
  c = mxGetPr(prhs[1]);
  c--;
  index = mxGetPr(prhs[2]);
  index--;
       
  plhs[0]=mxCreateDoubleMatrix(p,1,mxREAL);
  x=mxGetPr(plhs[0]);
  
  getx(x,u,c,index,m,p); 
}


