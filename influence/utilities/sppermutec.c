#include "mex.h"
#include <math.h>
// sppermutec Permutes the ordering of the variables in a virtually
// multidimensional sparse array.
// The logic of this algorithm is as follows:
// get the output row and column indices of each non-zero element
// and count the number of elements in each column. Then create the
// output column JC vector. Then put the data and the row indices into
// the memory area associated with the associated column. Finally
// sort the data in each column so it is ordered by row.

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

#define STACKSIZE 64  // will never need more than this for 64 bit integers
#define SWITCH    16  // use insertionsort for segments less than SWITCH

void insertionsort(mwSize *a, double *v, mwSize n){
  mwSize i, j, ax;
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
void quicksort(mwSize *a, double *v, mwSize n) {
  mwSize  beg[STACKSIZE], end[STACKSIZE], L, R, b0, e0, sptr=0,  m;
  mwSize apiv, atemp;
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
  mwSize *Irx, *Irxptr, *Irxend, *Jcx, *Ry, *Cy, *order, sep, d;
  double *x, *y;
  mwSize *Iry, *Jcy, *count, *Jcyi, *nx, *xfact, *yrfact, *ycfact;
  mwSize  rout, cout, ry, cy, ix, cx, rx, ryj, cyj, temp, temp0, temp1;
  mwSize i, j, nnz, k, colsx, rowsx, oi, oi1, kend;

  if (nrhs<4 || nrhs>4)
     mexErrMsgTxt("Four parameters must be passed");
  if (nlhs>1)
     mexErrMsgTxt("Only one output is created");
  for (j=0; j<4; j++){
    if (!mxIsDouble(prhs[j]))
      mexErrMsgTxt("Input must all be double arrays"); 
    if (mxIsComplex(prhs[j]))
      mexErrMsgTxt("Inputs cannot be complex");
    if (j>0 && mxIsSparse(prhs[j]))
      mexErrMsgTxt("Inputs 2-4 must be full");
  }
  if (!mxIsSparse(prhs[0]))
      mexErrMsgTxt("First input must be sparse");

d=(mwSize)mxGetNumberOfElements(prhs[1]);
 
if (d!=mxGetNumberOfElements(prhs[2]))
  mexErrMsgTxt("order and size vectors must be the same length");
if (mxGetNumberOfElements(prhs[3])!=2)
  mexErrMsgTxt("output size must be a 2-vector");

rowsx=(mwSize)mxGetM(prhs[0]);
colsx=(mwSize)mxGetN(prhs[0]);
x  = mxGetPr(prhs[0]);
Irx=mxGetIr(prhs[0]);
Jcx=mxGetJc(prhs[0]); 
nnz=(mwSize)Jcx[colsx];
rout=(mwSize) mxGetPr(prhs[3])[0];
cout=(mwSize) mxGetPr(prhs[3])[1];

nx=mxCalloc(6*d,sizeof(mwSize));
order=nx+d;
xfact=order+d;
yrfact=xfact+d;
ycfact=yrfact+d;

for (i=0;i<d;i++){
  order[i]=(mwSize)(mxGetPr(prhs[1])[i])-1;  // convert to 0-base
     nx[i]=(mwSize)(mxGetPr(prhs[2])[i]); 
 yrfact[i]=0;
}

plhs[0]=mxCreateSparse(rout, cout, nnz, mxREAL);
y  =mxGetPr(plhs[0]);
Iry=mxGetIr(plhs[0]);
Jcy=mxGetJc(plhs[0]); 
count=Jcy+1;

oi1=order[0];
if (rout>1){
  sep=0;
  yrfact[oi1]=1; 
}
else{
  sep=2;
  ycfact[oi1]=1; 
}
for (i=0;i<d-1;){
  oi=oi1;
  oi1=order[++i];
  if (sep==0){
    temp=yrfact[oi]*nx[oi];
    if (temp==rout){
      sep=1;
      ycfact[oi1]=1;
    }
    else{
      yrfact[oi1]=temp;
    }
  }
  else
    ycfact[oi1]=ycfact[oi]*nx[oi];
}

xfact[0]=1; sep=0;
for (i=0;i<d-1;){
  temp=xfact[i]*nx[i];
  i++;
  if (sep==0 && temp==rowsx){
    sep=i;
    xfact[i]=1;
  }
  else xfact[i]=temp;
}

if (sep==0){
  if (colsx==1){ sep=d;}
  else { if (rowsx>1){
    mexErrMsgTxt("dimensions of x must be aligned to rows and columns of x");}}
}

// row and column indices for the output
Ry=mxCalloc(2*nnz,sizeof(mwSize));
Cy=Ry+nnz;

//for (i=0;i<d;i++)printf("%3i  %3i  %3i\n",xfact[i],yrfact[i],ycfact[i]);

Irxend=Irx;
Irxptr=Irx;
k=0;
kend=0;
for (j=0;j<colsx;){
  cx=j;
  if (kend<Jcx[++j]){  //loop over columns
    // get output row and column indices associated with first element of column j
    kend=Jcx[j];
    ryj=0; cyj=0;
    if (sep<d){
      for (i=d-1;i>sep;i--){
        ix=cx/xfact[i];
        cx-=ix*xfact[i];
        if (yrfact[i]>0) ryj += ix*yrfact[i];
        else             cyj += ix*ycfact[i];
      }
      if (yrfact[sep]>0) ryj += cx*yrfact[sep];
      else               cyj += cx*ycfact[sep];
    }
    // loop over non-zeros in column j
    for (;k<kend;){
      ry=ryj; cy=cyj;
      if (sep>0){
        rx=Irx[k];
        for (i=sep-1;i>0;i--){
          ix=rx/xfact[i];
          rx-=ix*xfact[i];
          if (yrfact[i]>0) ry += ix*yrfact[i];
          else             cy += ix*ycfact[i];
        }
        if (yrfact[0]>0) ry += rx*yrfact[0];
        else             cy += rx*ycfact[0];
      }
      // store the row and column indices of element k and 
      // increment the column counter
      /* check on size of indices
      if (ry>=rout){
         printf("%1i %1i\n",ry,rout);
         mexErrMsgTxt("row index too big");
      }
      if (cy>=cout){
         printf("%1i %1i\n",cy,cout);
         mexErrMsgTxt("column index too big");
      }
      */
      Ry[k]=ry;
      Cy[k++]=cy;
      count[cy]++;
    }
  }
}

// get the cumulative column counts
temp1=0;
for (i=1;i<=cout;i++) {
  temp0=temp1;
  temp1=Jcy[i];
  Jcy[i]=Jcy[i-1]+temp0;
}

// put the values and row indices into the area of memory
// set aside for appropriate column
Jcy++;
for (i=0;i<nnz;i++) {
  Jcyi=Jcy+Cy[i];
  Iry[*Jcyi]=Ry[i];
  y[(*Jcyi)++]=x[i];
}
Jcy--;

// sort the data in each column
for (i=0;i<cout;) {
  k=Jcy[i]; j=Jcy[++i]-k;
  if (j<SWITCH)    insertionsort(Iry+k,y+k,j);
  else                 quicksort(Iry+k,y+k,j);
}

mxFree(Ry);
mxFree(nx);
}