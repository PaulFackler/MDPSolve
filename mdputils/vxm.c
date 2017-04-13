#include "mex.h"
#include <math.h>
/*
% devecxmat Computes diag(a)*B
% USAGE
%   C=diagmult(a,B);
% INPUTS
%   a : full m-vector (or scalar)
%   B : mxn matrix (full or sparse)
% OUTPUT
%   C   : mxn matrix 
%
% Note: not implemented for complex matrices or matrices 
% with data type other than double. a must be full but
% B can be sparse or full.

% Copyright (c) 2010, Paul L. Fackler, NCSU
% paul_fackler@ncsu.edu
*/


/* diag(a) times B - full B */
void daxbf(double *a, double *B, mwSize m, mwSize n)
{
  double *aptr, *aend, *Bend;
  aend=a+m;
  Bend=B+m*n;
  while (B<Bend){
    aptr=a;
    while (aptr<aend) *B++ *= *aptr++;
  }
}

/* diag(a) times B - sparse B */
void daxbs(double *a, double *B, mwIndex *Bi, mwSize nnz)
{
  double *Bend;
  Bend=B+nnz;
  while (B<Bend) *B++ *= a[*Bi++];
}




void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{  double *A, *B, *Bend, a;
   mwSize m, n, na;

   if (nrhs!=2)
      mexErrMsgTxt("Two parameters must be passed");
   if (nlhs>1)
      mexErrMsgTxt("Only one output is created");
   if (!mxIsDouble(prhs[0]) ||  mxIsSparse(prhs[0]))
      mexErrMsgTxt("First input must be a full vector"); 
   if (!mxIsDouble(prhs[1]) &&  !mxIsSparse(prhs[1]))
      mexErrMsgTxt("Second input must be numeric");
   if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
      mexErrMsgTxt("Inputs must be real");

   m=mxGetM(prhs[1]); 
   n=mxGetN(prhs[1]);
   na=mxGetNumberOfElements(prhs[0]);
   if (na==1){  /* a is scalar */
       a=*mxGetPr(prhs[0]);
       plhs[0]=mxDuplicateArray(prhs[1]);
       B=mxGetPr(plhs[0]);
       if (mxIsSparse(plhs[0])){
         mwIndex *Bj;
         Bj=mxGetJc(plhs[0]);
         m=Bj[n];  /* number of non-zero elements */
       }
       else m *= n;
       Bend=B+m;
       while (B<Bend) *B++ *= a;   
   }
   else if (na==m){ /* a has the right number of elements */
     A=mxGetPr(prhs[0]);
     plhs[0]=mxDuplicateArray(prhs[1]);
     B=mxGetPr(plhs[0]);
     if (mxIsSparse(plhs[0])){ 
       mwIndex *Bi, *Bj;
       Bi=mxGetIr(plhs[0]);
       Bj=mxGetJc(plhs[0]);
       daxbs(A,B,Bi,Bj[n]);
     }
     else
       daxbf(A,B,m,n);
   }
   else
      mexErrMsgTxt("Inputs are not conformable");
}
