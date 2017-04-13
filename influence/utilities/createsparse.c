#include "mex.h"
#include <math.h>

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

void insertionsort(mwSize *a, double *v, size_t n){
  size_t i, j, ax;
  double vx; 

  for (i=1; i < n; i++) {
    ax = a[i]; vx=v[i];
    j = i;
    while (j>0 && a[j-1] > ax ) {
      a[j] = a[j-1]; v[j] = v[j-1]; 
      j--;
    }
    a[j] = ax; v[j] = vx;
  }
}

// This version of quicksort works for the case when all the values of the sort
// key are distinct - it has not been tested for situations with duplicates.
// It sorts on a but rearranges v simultaneously.
void quicksort(size_t *a, double *v, size_t n) {
  #define STACKSIZE 64 // will never need more than this for 64 bit integers
  #define SWITCH 16    // use insertionsort for segments less than SWITCH
  size_t  beg[STACKSIZE], end[STACKSIZE], L, R, b0, e0, sptr=0,  m;
  size_t apiv, atemp;
  double vpiv, vtemp;
  beg[0]=0; end[0]=n-1;
  while (true) {
    b0=beg[sptr]; e0=end[sptr];
    L=b0; R=e0;
    if (R-L<SWITCH){
      insertionsort(a+L,v+L,R-L+1); // R-L is one less than the number of elements passed
      if (sptr==0) break;  // done
      sptr--;
    }
    else { 
      // get pivot using median of three and make appropriate swaps
      // the swapping puts the pivot in place 0, the value less than the
      // pivot in place 1 and the value greater than the pivot in place
      // n-1. The value in place 1 goes to place m. Then only places
      // 2 through n-2 need be partitioned.
      m=L+(R-L)/2;   // midpoint
      if (a[L]<a[m])
        if (a[m]<a[R])   // a[L]<a[m]<a[R]
            { apiv=a[m]; a[m]=a[L+1]; a[L+1]=a[L];
              vpiv=v[m]; v[m]=v[L+1]; v[L+1]=v[L];}
        else
          if (a[L]<a[R]) // a[L]<a[R]<a[m]
            { apiv=a[R]; a[R]=a[m]; a[m]=a[L+1]; a[L+1]=a[L];
              vpiv=v[R]; v[R]=v[m]; v[m]=v[L+1]; v[L+1]=v[L];}
          else           // a[R]<a[L]<a[m]
            { apiv=a[L]; a[L]=a[R]; a[R]=a[m]; a[m]=a[L+1]; a[L+1]=a[L];
              vpiv=v[L]; v[L]=v[R]; v[R]=v[m]; v[m]=v[L+1]; v[L+1]=v[L];}
      else
        if (a[L]<a[R])   // a[m]<a[L]<a[R]
            { apiv=a[L]; atemp=a[L+1]; a[L+1]=a[m]; a[m]=atemp;
              vpiv=v[L]; vtemp=v[L+1]; v[L+1]=v[m]; v[m]=vtemp;}
        else                
          if (a[m]<a[R]) // a[m]<a[R]<a[L]
            { apiv=a[R]; atemp=a[L+1]; a[L+1]=a[m]; a[m]=atemp; a[R]=a[L]; 
              vpiv=v[R]; vtemp=v[L+1]; v[L+1]=v[m]; v[m]=vtemp; v[R]=v[L];}
          else           // a[R]<a[m]<a[L]              
            { apiv=a[m]; a[m]=a[L+1]; a[L+1]=a[R]; a[R]=a[L];
              vpiv=v[m]; v[m]=v[L+1]; v[L+1]=v[R]; v[R]=v[L];} 
      L+=2; R--;  // reduce the partition space by 3 elements
      // partition
      while (true) {
        while (a[R] > apiv) R--; 
        while (a[L] < apiv) L++; 
        if (L>=R) break;
        atemp=a[L]; a[L]=a[R]; a[R]=atemp;
        vtemp=v[L]; v[L]=v[R]; v[R]=vtemp;
      }
      // swap the pivot to the rightmost left hand partition element
      a[b0]=a[R]; a[R]=apiv; 
      v[b0]=v[R]; v[R]=vpiv;
      R--;
      // now R is the rightmost element on the left and 
      //     L is the leftmost  element on the right
      // now sort the partitions
      if ((R-b0)>(e0-L)){
         // left partition is bigger - sort right first
        end[sptr++]=R; beg[sptr]=L; end[sptr]=e0;
      }
      else {
        // right partition is bigger - sort left first
        beg[sptr++]=L; beg[sptr]=b0; end[sptr]=R;
      }
    }
  }
}



void mexFunction(
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{ 
  double *Ir,*Jc, *Jcend;
  double *v, *y;
  mwSize *Iry, *Jcy, *count, *Jcyi, *iptr;
  mwSize  rout, cout;
  mwSize i, j, nnz, k;

  if (nrhs<5 || nrhs>5)
     mexErrMsgTxt("Five parameters must be passed");
  if (nlhs>1)
     mexErrMsgTxt("Only one output is created");
  for (j=0; j<5; j++){
    if (!mxIsDouble(prhs[j]))
      mexErrMsgTxt("Input must all be double arrays"); 
    if (mxIsComplex(prhs[j]))
      mexErrMsgTxt("Inputs cannot be complex");
    if (mxIsSparse(prhs[j]))
      mexErrMsgTxt("Inputs must be full");
  }

// get sizes of x 
nnz=mxGetNumberOfElements(prhs[0]);
if (mxGetNumberOfElements(prhs[1])!=nnz)
  mexErrMsgTxt("ry and cy must be the same length");
if (mxGetNumberOfElements(prhs[2])!=nnz)
  mexErrMsgTxt("ry and vy must be the same length");

Ir=mxGetData(prhs[0]);
Jc=mxGetData(prhs[1]);
v=mxGetPr(prhs[2]); 
rout=(mwSize)*mxGetPr(prhs[3]);
cout=(mwSize)*mxGetPr(prhs[4]);

if ((sizeof(mwSize)==32 && rout>2147483648/cout) || (sizeof(mwSize)==64 && rout>9007199254740992/cout))
  mexErrMsgTxt("matrix is too large to handle with linear indexing");

plhs[0]=mxCreateSparse(rout, cout, nnz, mxREAL);
y  =mxGetPr(plhs[0]);
Iry=mxGetIr(plhs[0]);
Jcy=mxGetJc(plhs[0]); 

count=mxCalloc(cout,sizeof(mwSize));

count--;
Jcend=Jc+nnz;
for (;Jc<Jcend;) {count[(mwSize)*Jc++]++;} Jc-=nnz;
for (i=1;i<cout;i++) {Jcy[i+1]=Jcy[i]+count[i];}
count++;
mxFree(count);
for (;Jc<Jcend;) {
  Jcyi=Jcy + (mwSize)*Jc++;
  Iry[*Jcyi]=(mwSize)*Ir++ - 1;
  y[(*Jcyi)++]=*v++;
}
for (iptr=Jcy+cout;Jcy<iptr;) {
  k=*Jcy++; j=*Jcy-k;
  if (j<SWITCH)    insertionsort(Iry+k,y+k,j);
  else                 quicksort(Iry+k,y+k,j);
}
}