#include "mex.h"
#include <math.h>
# include <time.h>
/*fs2f Computes products of tensors x and y with x and z full and y sparse
% USAGE
%   z=fs2f(X,x2z,Y,y2z);
% INPUTS
%   X     : a dx-dimensional tensor (implicitly)
%   x2z   : 1 or 2 row matrix with dx columns
%             x2z(1,j): the dimension of z associated with the jth dimension of X.
%               Use negative values to indicate that the dimension  should be summed
%             x2z(2,j): the size (number of elements) in the jth dimension of X
%               (if omitted it is assumed that the actual and implicit sizes are equal)
%   Y     : a dy-dimensional tensor (implicitly)
%   y2z   : 1 or 2 row matrix with dy columns (analogous to x2z)
% OUTPUT
%   z     : the resulting tensor
%
% Tensor multiplication takes three forms:
%   1) expansion - dimension multiplies all other dimensions
%   2) matching  - dimension multiplies only matching elelments in other array
%   3) summation - matches and sums over this dimension
% If both x2z and y2z contain a common element, this dimension is matched
% If the common element is negative it is also summed (which causes this
%   dimension of z to be a singleton and it can be squeezed out).
% Any value of x2z and y2z that is not in the other is an expansion dimension.
%
% For example kron(x,y) is obtained using x2z=[2 4] and y2z=[1 3].
% Ordinary matrix multiplication XY sets x2z=[1 -3] and y2z=[-3 2].
% Column-wise Kronecker products use x2z=[1 3] and y2z=[2 3].
% In each case nzout can be set if a matrix output is desired.
%
% It is not necessary for a tensor to be multi-dimensional (for sparse arrays
%  only 2-D are supported in MATLAB). Implicit dimensions can be used instead.
% For example, suppose that X is 20 x 15. It could implicitly be a 4-D array 
%   with implicit sizes nx=[5 4 3 5]; it must be the case that prod(nx)=numel(X).
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

#define abs(x) ((x)>=0 ? (x) : -(x))

// Ordered by y (the sparse matrix) and only non-zero elements are accessed. 
mxArray *fs2fsum(
  double  *valx, double  *valy, 
  mwSize *Iry, mwSize *Jcy,
  mwSize  rowy, mwSize coly,
  int dx, int dy,
  mwSize *nx, mwSize *ny,
  int *x2z, int *y2z){ 
  mxArray *z;
  double *valz, valkx, valky;
  mwSize *cumx, *cumy, *cumz, *subx, *nux;
  mwSize *fpyx, *fpyz, *fpxx, *fmxx, *fpxz, *fmxz, *nz;
  mwSize  nux0, fpxx0, fpxz0;
  mwSize indy;
  double *ptrx, *ptrz;
  mwIndex i, k, stepy, jy, subyi, ky, kyend;
  int ix, iy, iz, dz, dux, dux1;
  bool cont;  
    
  dz=0;
  for (ix=0; ix<dx; ix++) if (abs(x2z[ix])>dz) dz=abs(x2z[ix]);
  for (iy=0; iy<dy; iy++) if (abs(y2z[iy])>dz) dz=abs(y2z[iy]);
  // get array memory 
  nz=mxCalloc(12*dz,sizeof(mwSize));
  subx=nz+dz;
  nux=subx+dz;
  cumx=nux+dz; 
  cumy=cumx+dz; 
  cumz=cumy+dz; 
  fpyx=cumz+dz; 
  fpyz=fpyx+dz;
  fpxx=fpyz+dz;
  fpxz=fpxx+dz;
  fmxx=fpxz+dz;
  fmxz=fmxx+dz;
  
  for (iz=0; iz<dz; iz++) nz[iz]=0;
  for (ix=0; ix<dx; ix++){
    k=abs(x2z[ix]);
    if (k==0)
      mexErrMsgTxt("0-dimension is undefined");
    k--; // convert to 0-indexing
    nz[k]=nx[ix];
    if (x2z[ix]<0) nux[k]=1;  // use temporarily to designate summed dimensions
  }
  for (iy=0; iy<dy; iy++) {
    k=abs(y2z[iy]);
    if (k==0)
      mexErrMsgTxt("0-dimension is undefined");
    k--;  // convert to 0-indexing
    if (nz[k]>0 && nz[k]!=ny[iy]){
      mexErrMsgTxt("sizes of matched dimensions do not match");
    }
    if (nux[k]==1 && y2z[iy]>=0)
      mexErrMsgTxt("summed dimensions do not match");
    nz[k]=ny[iy];
  }
  // set unallocated and summed dimension sizes to 1
  for (iz=0; iz<dz; iz++) if (nz[iz]==0 || nux[iz]==1) nz[iz]=1;
  // construct control vectors
  cumy[0]=1; for (iy=1;iy<dy;iy++) cumy[iy]=cumy[iy-1]*ny[iy-1];
  cumx[0]=1; for (ix=1;ix<dx;ix++) cumx[ix]=cumx[ix-1]*nx[ix-1];
  cumz[0]=1; for (iz=1;iz<dz;iz++) cumz[iz]=cumz[iz-1]*nz[iz-1];
  // get controls for dimensions of y
  for (iy=0; iy<dy; iy++) {
    if (y2z[iy]>0) fpyz[iy]=cumz[y2z[iy]-1];
    else           fpyz[iy]=0;
    fpyx[iy]=0;
    for (ix=0;ix<dx; ix++){
      if (x2z[ix]==y2z[iy]) {fpyx[iy]=cumx[ix]; break;}
    }
  }
  // get controls for unique dimensions of x
  dux=0; 
  for (ix=0; ix<dx; ix++) {
    cont=true;
    // determine if ix is a unique dimension
    for (iy=0; iy<dy; iy++) if (x2z[ix]==y2z[iy]) {cont=false; break;}
    if (cont){
      nux[dux]=nx[ix];
      fpxx[dux]=cumx[ix];
      fmxx[dux]=fpxx[dux]*(nx[ix]-1);
      fpxz[dux]=cumz[x2z[ix]-1];
      fmxz[dux]=fpxz[dux]*(nz[x2z[ix]-1]-1);
      dux++;
    }
  }
  dux1=dux-1;
  // use when x has a single unique dimension
  if (dux==1) {nux0=nux[0]; fpxx0=fpxx[0]; fpxz0=fpxz[0];}
  z=mxCreateNumericArray(dz,nz,mxDOUBLE_CLASS, mxREAL);
  valz=mxGetPr(z);

  stepy=0; 
  // loop over elements of y
  for (jy=0; jy<coly; jy++, stepy+=rowy){  
    for (ky=Jcy[jy],kyend=Jcy[jy+1]; ky<kyend; ky++){
      indy=Iry[ky]+stepy;
      valky=valy[ky];
      // determine the subscripts for y associated with this element
      // and adjust the pointers to x and z accordingly
      ptrx=valx; ptrz=valz;
      for (iy=dy;iy>1;){
        iy--;
        subyi = indy/cumy[iy];  
        indy -= subyi*cumy[iy];
        ptrx += subyi*fpyx[iy];
        ptrz += subyi*fpyz[iy];
      }
      ptrx += indy * *fpyx;
      ptrz += indy * *fpyz;
      // loop over unique dimensions of x
      if (dux==0){
        *ptrz += *ptrx*valky;
      }
      else if (dux==1){  // special case with only a single unique dimension
        for (i=0; i<nux0; i++){
          *ptrz += *ptrx*valky;
          ptrx += fpxx0;
          ptrz += fpxz0;
        }
      }  // number of unique dimensions in x > 1
      else{
        for (ix=0; ix<dux; ix++) subx[ix]=0;
        cont=true;
        while (cont) {
          *ptrz += *ptrx*valky;
          for (ix=0; ix<dux; ix++){
            subx[ix]++;
            if (subx[ix]>=nux[ix]){
              if (ix>=dux1) cont=false; 
              else{
                subx[ix]=0;
                ptrx-=fmxx[ix];
                ptrz-=fmxz[ix];
              }
            }
            else{
              ptrx+=fpxx[ix];
              ptrz+=fpxz[ix];
              break;
            }
          }
        }
      }
    }
  }
  mxFree(nz);
  return(z); 
}

// Ordered by y (the sparse matrix) and only non-zero elements are accessed. 
mxArray *fs2f(
  double  *valx, double  *valy, 
  mwSize *Iry, mwSize *Jcy,
  mwSize  rowy, mwSize coly,
  int dx, int dy,
  mwSize *nx, mwSize *ny,
  int *x2z, int *y2z){ 
  mxArray *z;
  double *valz, valkx, valky;
  mwSize *cumx, *cumy, *cumz, *subx, *nux;
  mwSize *fpyx, *fpyz, *fpxx, *fmxx, *fpxz, *fmxz, *nz;
  mwSize indx, indy, indz;
  mwIndex i, k, stepy, jy, subyi, ky, kyend;
  int ix, iy, iz, dz, dux;
  bool cont;  
  dz=0;
  for (ix=0; ix<dx; ix++) if (abs(x2z[ix])>dz) dz=abs(x2z[ix]);
  for (iy=0; iy<dy; iy++) if (abs(y2z[iy])>dz) dz=abs(y2z[iy]);
  // get array memory 
  nz=mxCalloc(12*dz,sizeof(mwSize));
  subx=nz+dz;
  nux=subx+dz;
  cumx=nux+dz; 
  cumy=cumx+dz; 
  cumz=cumy+dz; 
  fpyx=cumz+dz; 
  fpyz=fpyx+dz;
  fpxx=fpyz+dz;
  fpxz=fpxx+dz;
  fmxx=fpxz+dz;
  fmxz=fmxx+dz;
  
  for (iz=0; iz<dz; iz++) nz[iz]=0;
  for (ix=0; ix<dx; ix++){
    k=abs(x2z[ix]);
    if (k==0)
      mexErrMsgTxt("0-dimension is undefined");
    k--; // convert to 0-indexing
    nz[k]=nx[ix];
  }
  for (iy=0; iy<dy; iy++) {
    k=abs(y2z[iy]);
    if (k==0)
      mexErrMsgTxt("0-dimension is undefined");
    k--;  // convert to 0-indexing
    if (nz[k]>0 && nz[k]!=ny[iy]){
      mexErrMsgTxt("sizes of matched dimensions do not match");
    }
    nz[k]=ny[iy];
  }
  // set unallocated and summed dimension sizes to 1
  for (iz=0; iz<dz; iz++) if (nz[iz]==0) nz[iz]=1;
  // construct control vectors
  cumy[0]=1; for (iy=1;iy<dy;iy++) cumy[iy]=cumy[iy-1]*ny[iy-1];
  cumx[0]=1; for (ix=1;ix<dx;ix++) cumx[ix]=cumx[ix-1]*nx[ix-1];
  cumz[0]=1; for (iz=1;iz<dz;iz++) cumz[iz]=cumz[iz-1]*nz[iz-1];
  // get controls for dimensions of y
  for (iy=0; iy<dy; iy++) {
    if (y2z[iy]>0) fpyz[iy]=cumz[y2z[iy]-1];
    else           fpyz[iy]=0;
    fpyx[iy]=0;
    for (ix=0;ix<dx; ix++){
      if (x2z[ix]==y2z[iy]) {fpyx[iy]=cumx[ix]; break;}
    }
  }
  // get controls for unique dimensions of x
  dux=0; 
  for (ix=0; ix<dx; ix++) {
    cont=true;
    // determine if ix is a unique dimension
    for (iy=0; iy<dy; iy++) if (x2z[ix]==y2z[iy]) {cont=false; break;}
    if (cont){
      nux[dux]=nx[ix];
      fpxx[dux]=cumx[ix];
      fpxz[dux]=cumz[x2z[ix]-1];
      fmxx[dux]=fpxx[dux]*(nx[ix]-1);
      fmxz[dux]=fpxz[dux]*(nz[x2z[ix]-1]-1);
      dux++;
    }
  }
  z=mxCreateNumericArray(dz,nz,mxDOUBLE_CLASS, mxREAL);
  valz=mxGetPr(z);
  stepy=0; 
  // loop over elements of y
  for (jy=0; jy<coly; jy++, stepy+=rowy){  
    for (ky=Jcy[jy],kyend=Jcy[jy+1]; ky<kyend; ky++){
      indy=Iry[ky]+stepy;
      valky=valy[ky];
      // determine the subscripts for z associated with this element
      indx=0; indz=0;
      for (iy=dy;iy>1;){
        iy--;
        subyi=indy/cumy[iy];  // suby[i]
        indx+=subyi*fpyx[iy];
        indz+=subyi*fpyz[iy];
        indy-=subyi*cumy[iy];
      }
      indx+=indy*fpyx[0];
      indz+=indy*fpyz[0];
      // loop over unique dimensions of x
      if (dux==0)
        valz[indz]=valx[indx]*valky;
      else{
        for (ix=0; ix<dux; ix++) subx[ix]=0;
        cont=true;
        while (cont) {
          valz[indz]=valx[indx]*valky;
          for (ix=0; ix<dux; ix++){
            subx[ix]++;
            if (subx[ix]>=nux[ix]){
              if ((ix+1)>=dux) {
                cont=false; 
              }
              else{
                subx[ix]=0;
                indx-=fmxx[ix];
                indz-=fmxz[ix];
              }
            }
            else{
              indx+=fpxx[ix];
              indz+=fpxz[ix];
              break;
            }
          }
        }
      }
    }
  }
  mxFree(nz);
  return(z); 
}


void mexFunction(
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{ 
  double *x, *y, *xind, *yind;
  mwSize *Iry, *Jcy;
  mwSize  rowy, coly;
  mwSize *nx, *ny, *nz;
  int dx, dy, *x2z, *y2z;
  int i,j;
  bool aresums;

  if (nrhs<4 || nrhs>4)
     mexErrMsgTxt("Four parameters must be passed");
  if (nlhs>1)
     mexErrMsgTxt("Only one output is created");
  for (j=0; j<4; j++){
    if (!mxIsDouble(prhs[j]))
      mexErrMsgTxt("Input must all be double arrays"); 
    if (mxIsComplex(prhs[j]))
      mexErrMsgTxt("Inputs cannot be complex");
  }
  if (mxIsSparse(prhs[0]))
      mexErrMsgTxt("x must be full (dense)");
  if (!mxIsSparse(prhs[2]))
      mexErrMsgTxt("y must be sparse");

// get dimensions of x and y
if (mxGetM(prhs[1])==2)  dx=(int)mxGetN(prhs[1]);
else                     dx=(int)mxGetNumberOfElements(prhs[1]);
if (mxGetM(prhs[3])==2)  dy=(int)mxGetN(prhs[3]);
else                     dy=(int)mxGetNumberOfElements(prhs[3]);

x2z=mxCalloc(dx+dy,sizeof(int));
y2z=x2z+dx;
// nz is a pointer placeholder because nx and ny might be redefined
nz=mxCalloc(1+dx+dy,sizeof(mwSize));
nx=nz+1;
ny=nx+dx;

xind=mxGetPr(prhs[1]);
yind=mxGetPr(prhs[3]);
aresums=false;
// get control information for x
if (mxGetM(prhs[1])==2){
  for (i=0;i<dx;i++){
    x2z[i]=(int)xind[2*i];
    if (x2z[i]<0) aresums=true;
    nx[i]=(mwSize)xind[2*i+1];
  }
}
else{
  if (mxGetNumberOfDimensions(prhs[0])!=dx)
    mexErrMsgTxt("Dimensions of X do not match x2z");
  nx=mxGetDimensions(prhs[0]);
  for (i=0;i<dx;i++){
    x2z[i]=(int)xind[i];
    if (x2z[i]<0) aresums=true;
  }   
}
// get control information for y
if (mxGetM(prhs[3])==2){
  for (i=0;i<dy;i++){
    y2z[i]=(int)yind[2*i];
    if (y2z[i]<0) aresums=true;
    ny[i]=(mwSize)yind[2*i+1];
  }
}
else{
  if (mxGetNumberOfDimensions(prhs[2])!=dy)
    mexErrMsgTxt("Dimensions of Y do not match y2z");
  ny=mxGetDimensions(prhs[2]);
  for (i=0;i<dy;i++){
    y2z[i]=(int)yind[i];
    if (y2z[i]<0) aresums=true;
  }   
}
 
x=mxGetPr(prhs[0]);
y=mxGetPr(prhs[2]);
Iry=mxGetIr(prhs[2]);
Jcy=mxGetJc(prhs[2]);
rowy=mxGetM(prhs[2]);
coly=mxGetN(prhs[2]);
if (aresums)
  plhs[0]=fs2fsum(x,y,Iry,Jcy,rowy,coly,dx,dy,nx,ny,x2z,y2z);
else
  plhs[0]=fs2f(x,y,Iry,Jcy,rowy,coly,dx,dy,nx,ny,x2z,y2z);
mxFree(nz);
mxFree(x2z);
}
