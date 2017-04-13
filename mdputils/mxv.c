#include "mex.h"
#include <math.h>
/*
% mxv Computes A*diag(b) (matrix times vector)
% USAGE
%   C=mxv(A,b);
% INPUTS
%   A : mxn matrix (full or sparse)
%   b : full n-vector (or scalar)
% OUTPUT
%   C   : mxn matrix 
%
% Note: not implemented for complex matrices or matrices 
% with data type other than double. b must be full but
% A can be sparse or full.

% Copyright (c) 2010, Paul L. Fackler, NCSU
% paul_fackler@ncsu.edu
*/


/* A times diag(b) - full A */
void axdbf(double *A, double *b, mwSize m, mwSize n)
{
  double *Ajend,*Aend, bval;
  Ajend=A;
  Aend=A+m*n;
  while(A<Aend){
    bval=*b++;
    Ajend += m;
    while(A<Ajend) *A++ *= bval;
  }
}

/* A times diag(b) - sparse A */
void axdbs(double *A, double *b, mwIndex *Aj, mwSize nnz)
{
  double *Aptr, *Anext, *Aend, bval;
  Aptr=A;
  Aend=A+nnz;
  while (Aptr<Aend){
    bval=*b++;
    Anext = A+ *(++Aj);
    while (Aptr<Anext) *Aptr++ *= bval;
  }
}



void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{  double *A, *B, b;
   mwSize m, n, nb;
   mwIndex i;
   bool overwrite;

   if (nrhs<2)
      mexErrMsgTxt("Two parameters must be passed");
   if (nrhs>3)
      mexErrMsgTxt("At most three parameters can be passed");
   if (nlhs>1)
      mexErrMsgTxt("Only one output is created");
   if (!mxIsDouble(prhs[0]) &&  !mxIsSparse(prhs[0]))
      mexErrMsgTxt("First input must be numeric"); 
   if (!mxIsDouble(prhs[1]) ||  mxIsSparse(prhs[1]))
      mexErrMsgTxt("Second input must be a full vector");
   if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
      mexErrMsgTxt("Inputs must be real");
   
   overwrite=false;
   if (nrhs>2 && *mxGetPr(prhs[2])!=0)  overwrite=true;
   
   m=mxGetM(prhs[0]); 
   n=mxGetN(prhs[0]);
   nb=mxGetNumberOfElements(prhs[1]);
   if (nb==1){  /* b is scalar */
       b=*mxGetPr(prhs[1]);
       if (overwrite) plhs[0]=prhs[0];
       else           plhs[0]=mxDuplicateArray(prhs[0]);
       A=mxGetPr(plhs[0]);
       if (mxIsSparse(plhs[0])){
         mwIndex *Aj;
         Aj=mxGetJc(plhs[0]);
         m=Aj[n];
       }
       else m *= n;
       for (i=0; i<m; i++) A[i]*=b;   
   }
   else if (nb==n){ /* b has the right number of elements */
     B=mxGetPr(prhs[1]);
     if (overwrite) plhs[0]=prhs[0];
     else           plhs[0]=mxDuplicateArray(prhs[0]);
     A=mxGetPr(plhs[0]);
     if (mxIsSparse(plhs[0])){
       mwIndex *Aj;
       Aj=mxGetJc(plhs[0]);
       axdbs(A,B,Aj,Aj[n]);
     }
     else
       axdbf(A,B,m,n);
   }
   else
      mexErrMsgTxt("Inputs are not conformable");
}
