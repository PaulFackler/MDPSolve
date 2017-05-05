#include "mex.h"
#include <math.h>
/*
dompurge purges dominated vectors - used for POMDPs
v=dompurge(v);
*/

mwSize purge(double *v, mwIndex *ind, mwSize n, mwSize m){
  mwIndex i, j, k, s;
  mwSize M;
  double d, tol, *vj, *vk;
  bool ambiguous;
  
  tol=5e-12;
  for (j=0; j<m; j++) ind[j]=j+1;
  M=m;
  vj=v;
  for (j=0; j<m-1; j++, vj+=n){
    if (ind[j]>0){
      vk=vj+n;
      for (k=j+1; k<m; k++, vk+=n){
        if (ind[k]>0){
          s=0; ambiguous=false;
          for (i=0; i<n; i++){
            d=vj[i] - vk[i];
            if (d>tol)
              if (s>=0) s++;
              else ambiguous=true;
            else if (d<-tol)
              if (s<=0) s--;
              else ambiguous=true;
            if (ambiguous) break;
          }
          if (!ambiguous)
            if (s<=0) {ind[j]=0; M--;}
            else      {ind[k]=0; M--;}
          if (ind[j]==0) break;
        }
      }
    }
  }
  return(M);
}
  
  

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double  *v, *Ind;
  mwSize  n, m, M;
  mwIndex i, j, *ind;

  /* Error checking on inputs */  
  if (nrhs!=1) mexErrMsgTxt("Incorrect number of input arguments");
  if (!mxIsDouble(prhs[0])) 
      mexErrMsgTxt("V must be double");
  if (mxIsSparse(prhs[0]))
      mexErrMsgTxt("V must be full (dense)");
  if (mxIsComplex(prhs[0]))
      mexErrMsgTxt("V must be real");
  
  n=mxGetM(prhs[0]);
  m=mxGetN(prhs[0]);
  v=mxGetPr(prhs[0]);
  
  ind=mxCalloc(m,sizeof(mwIndex));
  M=purge(v,ind,n,m);
  
  if (M<=0){
    printf("%2i \n",M); 
    for (i=0; i<m; i++) printf("%2i ",ind[i]); printf("\n");
    mexErrMsgTxt("M must be positive");
  }
        
  plhs[0]=mxCreateDoubleMatrix(1,M,mxREAL);
  Ind=mxGetPr(plhs[0]);
  for (j=0,i=0; i<m; i++) if (ind[i]>0) Ind[j++]=(double) i+1;
  mxFree(ind);
}


