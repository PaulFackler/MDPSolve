#include "mex.h"
#include <math.h>
/*
% spnormalize Normalizes so elements of specified dimensions sum to 1
% USAGE
%    P=spnormalize(P,n,sumvars);
% INPUTS
%    P  : intrinsically multidimensional array with d dimensions
%    n  : d-vector of dimension sizes (prod(n)=numel(P))
%    sumvars : list of dimensions on which to normalize
% OUTPUT
%    P  : normalized array of the same size as P
*/

/*
% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2014, Paul L. Fackler (paul_fackler@ncsu.edu)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without  
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%        this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice, 
%        this list of conditions and the following disclaimer in the 
%        documentation and/or other materials provided with the distribution.
%    * Neither the name of the North Carolina State University nor of Paul L. 
%        Fackler may be used to endorse or promote products derived from this 
%        software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% For more information, see the Open Source Initiative OSI site:
%   http://www.opensource.org/licenses/bsd-license.php
*/

extern mxArray *mxCreateSharedDataCopy(const mxArray *pr);

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double *P, *ps;
  mwSize *PIr, *PJc, *n, *cnp, *cns, *inds, *sumvar;
  mwSize nrow, ncol, j, jj, k, i, knext, d, q, nnz, indp, si, indsk, split, isj,maxi;
  bool overwrite;

   /* Error checking on inputs */  
  if (nrhs<3 || nrhs>4) 
      mexErrMsgTxt("Incorrect number of input arguments.");
  if (!mxIsDouble(prhs[0]))
      mexErrMsgTxt("Inputs must be double");
  if (!mxIsSparse(prhs[0]))
      mexErrMsgTxt("First input must be sparse");
  if (!mxIsDouble(prhs[1]))
      mexErrMsgTxt("Inputs must be double");
  if (!mxIsDouble(prhs[2]))
      mexErrMsgTxt("Inputs must be double");
  
  overwrite=false;
  if (nrhs==4)
    if (*mxGetPr(prhs[3])!=0) overwrite=true;

  d=(mwSize)mxGetNumberOfElements(prhs[1]);
  q=(mwSize)mxGetNumberOfElements(prhs[2]);
  if (q>d)
    mexErrMsgTxt("more summed variables than dimensions in P");
  
  P  =mxGetPr(prhs[0]);
  PIr=mxGetIr(prhs[0]);
  PJc=mxGetJc(prhs[0]);
  
  nrow=(mwSize)mxGetM(prhs[0]);
  ncol=(mwSize)mxGetN(prhs[0]);
  nnz =(mwSize)PJc[ncol];
  
  n=mxCalloc(4*d+nnz,sizeof(mwSize));
  sumvar=n+d;
  cnp=sumvar+d;
  cns=cnp+d;
  inds=cns+d;
  
  ps=mxGetPr(prhs[1]);
  for (i=0;i<d;i++) {
    n[i]=(mwSize)ps[i]; 
    sumvar[i]=0;
  }
  ps=mxGetPr(prhs[2]);
  for (i=0;i<q;i++) {
    jj=(mwSize)ps[i]-1;
    if (jj>=d)
      mexErrMsgTxt("illegal value for a summed variable");
    if (sumvar[jj]==1)
      mexErrMsgTxt("repeated values in sumvar are not allowed");
    sumvar[jj]=1;
  }
  
  cnp[0]=1;
  cns[0]=1;
  for (i=1;i<d;i++) {
    cnp[i]=cnp[i-1]*n[i-1]; 
    if (sumvar[i-1]) cns[i]=cns[i-1];
    else             cns[i]=cns[i-1]*n[i-1];
  }
  if (cnp[d-1]*n[d-1] > nrow*ncol)
      mexErrMsgTxt("prod(n) too big");
  
  // determine the maximum number of summed elements
  if (d>1){
    if (sumvar[d-1]==0) maxi=cns[d-1]*n[d-1];
    else                maxi=cns[d-1];
  }
  else{
    if (sumvar[0]==0) maxi=n[0];
    else              maxi=1;
  }
  // determine if split falls at a boundary
  split=d;
  for (i=0;i<d;i++) if (cnp[i]==nrow) {split=i; break;}
  if (split==d && ncol>1)
      mexErrMsgTxt("matrix not split correctly");
  for (i=split;i<d;i++) {cnp[i]/=nrow;}
  
  //  if requested set up second output, otherwise use tempory memory
  if (nlhs>1){
    plhs[1]=mxCreateDoubleMatrix(cns[split],maxi/cns[split],mxREAL);
    ps=mxGetPr(plhs[1]);
  }
  else
    ps=mxCalloc(maxi,sizeof(double));
   
  jj=0;
  k=0;
  for (j=0;j<ncol;j++,jj+=nrow){  // loop over columns
    knext=PJc[j+1];
    if (k<knext){                 // if column is empty then skip
      isj=0;
      if (split<d){
        indp=j;
        for (i=d-1;i>split;i--){
          si=indp/cnp[i];
          indp -= si*cnp[i];
          if (sumvar[i]==0) isj += si*cns[i];
        }
        if (sumvar[split]==0) isj += indp*cns[split];
      }
      // loop over non-zero elements in the jth column
      for (;k<knext;k++){
        indsk=isj;
        if (split>0){
          indp=PIr[k];
          for (i=split-1;i>0;i--){
            si=indp/cnp[i];
            indp -= si*cnp[i];
            if (sumvar[i]==0) indsk += si*cns[i];
          }
          if (sumvar[0]==0) indsk += indp;
        }
        ps[indsk]+=P[k];  // sum into appropriate  element of ps
        inds[k]=indsk;    // store the linear index in ps
      }
    }
  }
  // plhs[0]=prhs[0];
  // uses an undocumented feature described here:
  // http://www.mk.tu-berlin.de/Members/Benjamin/mex_sharedArrays
  if (overwrite)
    plhs[0]=mxCreateSharedDataCopy(prhs[0]);
  else
    plhs[0]=mxDuplicateArray(prhs[0]);
  P=mxGetPr(plhs[0]);
  for (k=0; k<nnz; k++){
    P[k]/=ps[inds[k]];
  }
  if (nlhs==0) mxFree(ps);
  mxFree(n);
}