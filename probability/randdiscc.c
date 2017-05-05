#include "mex.h"
#include <math.h>

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double *p, *u, *ind, uj, *x, *pptr;
  mwSize q, m, n;
  mwIndex i, j, k, pind, pend, *Ir, *Jc;

  /* Error checking on inputs */  
  
  if (nrhs<2 || nrhs>3) 
      mexErrMsgTxt("Incorrect number of input arguments.");
  
  for (i=0; i<nrhs; i++){
    if (!mxIsDouble(prhs[i]))
        mexErrMsgTxt("Inputs must be double");
  }
  if (mxIsSparse(prhs[1]))
        mexErrMsgTxt("random variates must be dense");
             
  p=mxGetPr(prhs[0]);
  m=(mwSize)mxGetM(prhs[0]);
  n=(mwSize)mxGetN(prhs[0]);
  u=mxGetPr(prhs[1]);
  q=(mwSize)mxGetNumberOfElements(prhs[1]);
  if (nrhs==3){
    if (mxIsSparse(prhs[2]))
        mexErrMsgTxt("index vector must be dense");
    if (mxGetNumberOfElements(prhs[2])!=q)
        mexErrMsgTxt("u and index must have the same # of elements");
    ind=mxGetPr(prhs[2]);
  }
  else{
    if (n!=q)
        mexErrMsgTxt("CPT must have the same # of columns as elements in u if no index specified");
  }
  plhs[0]=mxCreateDoubleMatrix(q,1,mxREAL);
  x=mxGetPr(plhs[0]);
  // p is sparse
  if (mxIsSparse(prhs[0])){
    Ir=mxGetIr(prhs[0]);
    Jc=mxGetJc(prhs[0]);
    for (j=0; j<q;){
      uj=u[j];
      if (nrhs==2){
        pind=Jc[j];
        pend=Jc[++j];
      }
      else{
        i=(mwIndex)ind[j++];
        if (i>n || i<1) mexErrMsgTxt("index out of bounds");
        pind=Jc[i-1];
        pend=Jc[i];
      }
      for (;pind<pend;pind++){
        uj -= p[pind];
        if (uj<=0){*x++ = Ir[pind]+1; break;}
      }
    }
  }
  // p is full
  else{
    for (j=0; j<q; j++){
      uj=u[j];
      if (nrhs==2) i=j;
      else {
        if (ind[j]>n || ind[j]<1) mexErrMsgTxt("index out of bounds");
        i=(mwIndex)ind[j]-1; 
      }
      pptr=p+i*m;
      for (k=0; k<m; ){
        uj -= pptr[k++];
        if (uj<=0){*x++ = k; break;}
      }
    }
  }
}
