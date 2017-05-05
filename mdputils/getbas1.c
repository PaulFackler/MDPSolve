#include "mex.h"
#include <math.h>

/* 1-D linear interpolation on a regular grid
 * Called by rectbas.m */
  
#define min(x,y) ((x)<=(y) ? (x) : (y))
#define max(x,y) ((x)>=(y) ? (x) : (y))


mwIndex lookup(double xi, double *table, mwSize n) {
   mwSignedIndex j, jlo, jhi;
   mwSize        inc, n1, n2;

  /* handle 1-value lists separately */
  if (n==1) {
    return(1);  
  }
  else {
    n1=n-1;
    n2=n-2;
    jlo=0;
    inc=1;
    if (xi>=table[jlo]) {
      jhi=jlo+1;
      while (xi>=table[jhi]) {
        jlo=jhi;
        jhi+=inc;
        if (jhi>=n) { jhi=n; break; }
        else { inc += inc; }
      }
    }
    else {
      jhi=jlo;
      jlo--;
      while (xi<table[jlo]) {
        jhi=jlo;
        jlo-=inc;
        if (jlo<0) { jlo=-1; break; }
        else { inc += inc; }
      }
    }
    while (jhi-jlo>1) {
      j=(jhi+jlo)/2;
      if (xi>=table[j]) jlo=j; 
      else jhi=j; 
    }
    return(jlo);
  }
}

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  mwSize N, n; 
  mwIndex i, j, *ir, *jc, indi;
  double *S, *s, *b, si, bi, factor;
  int ii;
  bool evenspacing;

  /* Error checking on inputs */  
  if (nrhs<3) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>3) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments.");
  for (ii=0; ii<nrhs; ii++)
  {
    if (!mxIsDouble(prhs[ii]) || mxIsSparse(prhs[ii]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[ii]))
      mexErrMsgTxt("X must be real.");
  }
  if (mxGetNumberOfElements(prhs[2])>1)
      mexErrMsgTxt("Third input must be scalar");
  
  s=mxGetPr(prhs[0]);
  S=mxGetPr(prhs[1]);
  n=mxGetNumberOfElements(prhs[0]);
  
  // maximal sizes of table
  if (sizeof(mwIndex)==32){
    if (n > 2147483647)
      mexErrMsgTxt("Not supported for tables of size greater than 2147483647");
  }
  else{
    if (n > 2305843009213693951)
      mexErrMsgTxt("Not supported for tables of size greater than 2305843009213693951");
  }
  if (n<2)
      mexErrMsgTxt("Table must contain at least 2 elements");
  N=mxGetNumberOfElements(prhs[1]);
  
  if (*mxGetPr(prhs[2])==0)
    evenspacing=false;
  else
    evenspacing=true;
  
    // allocate memory for outputs
  if (nlhs<2){
    plhs[0]=mxCreateSparse(n, N, 2*N, mxREAL);
    b = mxGetPr(plhs[0]);
    ir = mxGetIr(plhs[0]);
    jc = mxGetJc(plhs[0]);
    for (j=0;j<=N;j++) jc[j]=j*2;
  }
  else{
    plhs[0]=mxCreateDoubleMatrix(2,N,mxREAL);
    b = mxGetPr(plhs[0]);
    if (sizeof(mwIndex)==4)
      plhs[1]=mxCreateNumericMatrix(2,N,mxUINT32_CLASS,mxREAL);
    else
      plhs[1]=mxCreateNumericMatrix(2,N,mxUINT64_CLASS,mxREAL);
    ir = mxGetData(plhs[1]);
  }
  
  // loop over the input data
  if (evenspacing){
    factor=(n-1)/(s[n-1]-s[0]);
    for (i=0;i<N;i++){
      if (S[i]<=s[0]) indi=0;
      else indi=min(floor((S[i]-s[0])*factor),n-2);
      *ir++=indi;
      si=s[indi];
      bi=(S[i]-si)/(s[++indi]-si);
      *ir++=indi;
      *b++=1-bi; *b++=bi; 
    }
  }
  else{
    for (i=0;i<N;i++){
      if (S[i]<=s[0]) indi=0;
      else indi=min(lookup(S[i],s,n),n-2);
      *ir++=indi;
      si=s[indi];
      bi=(S[i]-si)/(s[++indi]-si);
      *ir++=indi;
      *b++=1-bi; *b++=bi; 
    }
  }
}

