#include "mex.h"
#include <math.h>

/*
% CROSSORTHMOM Computes the complete set of Hermite polynomials up to order k
% USAGE
%   m=crossorthmom(y,k);
% INPUTS
%   y   : an n x d matrix
%   k   : a scalar positive integer
% OUTPUTS
%   m : an n x p row matrix where p=(k+d)!/(k!d!)
% Example: 
%   If y is n x 2 and k=3, m is
%   [H0 H1(y1) H1(y2) H2(y1) H1(y1)H1(y2) H2(y2) ...
%          H3(y1) H2(y1)H1(y2) H1(y1)H2(y2) H3(y2)]
% where Hj is the modified Hermite polynomial
%    Hj(x) = (1/sqrt(j))xHj-1(x)-sqrt((j-1)/j)Hj-2(x)
% with H0=1 and H1(x)=x
*/

void crossorthmom(double *y,double *m,mwSize n,mwSize d,mwSize k){
mwSize i, j, jj, kk, count, ind1, ind2;
mwSize *locs;
double factor1, factor2, *ptrm, *ptri, *ptr1, *ptr2;
  for (i=0; i<n;) m[i++]=1;
  if (k>0) {
    if (n>0 && d>0){
      memcpy(m+n,y,n*d*sizeof(double));
      locs=mxCalloc(d*(k+1),sizeof(mwSize));
      for (i=0;i<d;i++) {locs[i]=0; locs[i+d]=(i+1)*n;}
      count=(d+1)*n;
      for (j=2; j<=k; j++){
        factor1=1/sqrt(j);
        factor2=sqrt(((double)j-1)/j);
        for (i=0; i<d; i++){
          locs[i+j*d]=count;
          ptri=y+i*n;
          ptr1=m+locs[i+(j-1)*d];
          ptr2=m+locs[i+(j-2)*d];
          ptrm=m+count;
          for (kk=0;kk<n;kk++) 
            {ptrm[kk]=ptri[kk]*ptr1[kk]*factor1-ptr2[kk]*factor2;}
          count+=n;
          if (i<(d-1)){
            for (jj=j-1; jj>=1; jj--){
              ptri=m+locs[i+jj*d];
              ind1=locs[i+1+(j-jj)*d];
              ind2=locs[d-1+(j-jj)*d];
              for (;ind1<=ind2; count+=n)
                for (ptrm=m+count, kk=0; kk<n; kk++) ptrm[kk]=ptri[kk]*m[ind1++];
            }
          }
        }
      }
      mxFree(locs);
    }
  }
}

void mexFunction(
   mwSize nlhs, mxArray *plhs[],
   mwSize nrhs, const mxArray *prhs[])
{
  mwSize n, d, k, p, i, denom;
  mwSize *locs;
  double *y, *m;
  if (nrhs<2) mexErrMsgTxt("At least two inputs are required.");
  if (nrhs>2) mexErrMsgTxt("Only two inputs are allowed.");
  if (nlhs>1) mexErrMsgTxt("Only one output is allowed.");
  for (i=0;i<nrhs;i++){
    if (!mxIsDouble(prhs[i]) || mxIsSparse(prhs[i]))
      mexErrMsgTxt("Inputs must be full");
    if (mxIsComplex(prhs[i]))
      mexErrMsgTxt("Input must be real");
  }
  n=mxGetM(prhs[0]);
  d=mxGetN(prhs[0]);
  if (mxGetNumberOfElements(prhs[1])!=1)
    mexErrMsgTxt("k must be a scalar");
  k=(mwSize)*mxGetPr(prhs[1]);
  if (k<0) mexErrMsgTxt("k must be a non-negative integer");
  /* Determine the size of the result */
  p=1;denom=1;
  if (k>d) for (i=1;i<=d;i++) {p*=k+i; denom*=i;}
  else     for (i=1;i<=k;i++) {p*=d+i; denom*=i;}
  p/=denom;
  y=mxGetPr(prhs[0]);
  plhs[0]=mxCreateDoubleMatrix(n,p,mxREAL);
  m=mxGetPr(plhs[0]);
  crossorthmom(y,m,n,d,k);
}
