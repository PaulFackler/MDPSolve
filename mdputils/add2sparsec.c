#include "mex.h"
#include <math.h>
/*
  add2sparsec(A,M,startcol)
  Places an m-row sparse matrix M into a larger m-row sparse matrix A
    starting at column startcol+1
  It must be the case that cols(M)+startcol<=cols(A)
  
  it must also be true that there is enough memory allocated to A to hold the elements of M
  currently no checks are made so use of this function could cause Matlab to crash
  use add2sparse instead.
*/

extern mxArray *mxCreateSharedDataCopy(const mxArray *pr);

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double *A, *M;
  mwSize mA, nA, mM, nM;
  mwIndex *irA, *jcA, *irM, *jcM, jA,  i, j, startcol, k, knext;

  /* Error checking on inputs */  
  if (nrhs<3 || nrhs>3) 
      mexErrMsgTxt("Incorrect number of input arguments.");
  if (!mxIsDouble(prhs[0]))
      mexErrMsgTxt("Inputs must be double");
  if (!mxIsDouble(prhs[1]))
      mexErrMsgTxt("Inputs must be double");
  if (!mxIsSparse(prhs[0]))
      mexErrMsgTxt("Input 1 must be sparse");
  if (!mxIsDouble(prhs[2]))
      mexErrMsgTxt("Inputs must be double");
  if (mxGetNumberOfElements(prhs[2])!=1)
      mexErrMsgTxt("Third input must be scalar");
              
  startcol=*mxGetPr(prhs[2]);
        
  mA=mxGetM(prhs[0]);
  mM=mxGetM(prhs[1]);
  if (mA!=mM)
     mexErrMsgTxt("Matrices must have the same number of rows");

  nA=mxGetN(prhs[0]);
  nM=mxGetN(prhs[1]);
  if (nA-startcol<nM){
     printf("n(A):%1i  nM+startcol: %1i \n",nA,startcol+nM);
     mexErrMsgTxt("Not enough columns in A to hold M");
  }
  
  // plhs[0]=prhs[0];
  // uses an undocumented feature described here:
  // http://www.mk.tu-berlin.de/Members/Benjamin/mex_sharedArrays
  plhs[0]=mxCreateSharedDataCopy(prhs[0]);
  
  irA=mxGetIr(plhs[0]);
  jcA=mxGetJc(plhs[0]);
  
  A=mxGetPr(plhs[0]);
  M=mxGetPr(prhs[1]);
  
  jA=jcA[startcol];
  if (mxIsSparse(prhs[1])){
    irM=mxGetIr(prhs[1]);
    jcM=mxGetJc(prhs[1]);
    if (jcM[nM]+jA>mxGetNzmax(prhs[0])){
      printf("nnz(M)+nnz(A): %1i  nzmax(A): %1i\n",jcM[nM]+jA,mxGetNzmax(plhs[0]));
      mexErrMsgTxt("Not enough memory allocated to A to hold M");
    }
    k=0;
    for (j=1;j<=nM;j++){
      knext=jcM[j];
      while (k<knext){
        A[jA]=M[k];
        irA[jA++]=irM[k++];
      }
      jcA[startcol+j]=jA;
    }
  }
  else {
    k=0;
    for (j=1;j<=nM;j++){
      for (i=0;i<mA;i++,k++){
        if (M[k]!=0){
          A[jA]=M[k];
          irA[jA++]=i;
        }
      }
      jcA[startcol+j]=jA; 
    }
  }
  for (j=startcol+nM+1; j<=nA;) jcA[j++]=jA;
}
