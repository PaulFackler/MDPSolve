#include "mex.h"
#include <math.h>

#ifdef old32
#define mwSize int
#define mwIndex int
#endif

/* 
KRON MEX replacement for kron
*/

void kronffc(double *ReA, double *ImA, mwSize mA, mwSize nA, 
             double *ReB, double *ImB, mwSize mB, mwSize nB, 
             double *ReC, double *ImC)
{
  mwIndex iA, iB, jA, jB;
  double Avalr, Avali;
  bool Acomp,Bcomp;
    
  if (ImA==NULL) {Acomp=false; Avali=0;} else Acomp=true;
  if (ImB==NULL) Bcomp=false; else Bcomp=true;

  for (jA=0; jA<nA; jA++){
    for (jB=0; jB<nB; jB++)
    {
      for (iA=0; iA<mA; iA++)
      {
        Avalr=ReA[iA]; if (Acomp) Avali=ImA[iA];
        if (Bcomp)
          for (iB=0; iB<mB; iB++, ReC++,ImC++){
            *ReC=Avalr*ReB[iB]-Avali*ImB[iB];
            *ImC=Avali*ReB[iB]+Avalr*ImB[iB];
          }
        else
          for (iB=0; iB<mB; iB++, ReC++,ImC++){
            *ReC=Avalr*ReB[iB];
            *ImC=Avali*ReB[iB];
          }
      }
      ReB+=mB;
      if (Bcomp) ImB+=mB;
    }
    ReA+=mA;
    if (Acomp) ImA+=mA;
    ReB-=mB*nB;
    if (Bcomp) ImB-=mB*nB;
  }
}

void kronsfc(double *ReA, double *ImA,  mwIndex *cindA, mwIndex *rindA, mwSize mA, mwSize nA, 
             double *ReB, double *ImB,                                  mwSize mB, mwSize nB, 
             double *ReC, double *ImC,  mwIndex *cindC, mwIndex *rindC)
{
  mwIndex iA, iB, jA, jB, jBmB, ck, Arow, Aend, ii;
  double Avalr, Avali;
  bool Acomp,Bcomp;
    
  if (ImA==NULL) {Acomp=false; Avali=0;} else Acomp=true;
  if (ImB==NULL) Bcomp=false; else Bcomp=true;

  ck=0;
  for (jA=0; jA<nA; jA++,cindA++)
    for (jB=0, jBmB=0; jB<nB; jB++,jBmB+=mB,*(++cindC)=ck)
      for (iA=*cindA,Aend=*(cindA+1); iA<Aend; iA++){
        Avalr=ReA[iA]; if (Acomp) Avali=ImA[iA]; 
        if (Avalr!=0 || Avali!=0){
          Arow=rindA[iA]*mB;
          if (Bcomp)
            for (iB=0, ii=jBmB; iB<mB; iB++, ii++){
              if (ReB[ii]!=0 || ImB[ii]!=0){
                *ReC++=Avalr*ReB[ii]-Avali*ImB[ii];
                *ImC++=Avalr*ImB[ii]+Avali*ReB[ii];
                *rindC++ = Arow+iB;
                ck++;
              }
            }
          else 
            for (iB=0, ii=jBmB; iB<mB; iB++, ii++){
              if (ReB[ii]!=0){
                *ReC++=Avalr*ReB[ii];
                *ImC++=Avali*ReB[ii];
                *rindC++ = Arow+iB;
                ck++;
              }
            }
        }
      }
}

void kronfsc(double *ReA, double *ImA,                                  mwSize mA, mwSize nA,
             double *ReB, double *ImB,  mwIndex *cindB, mwIndex *rindB, mwSize mB, mwSize nB,
             double *ReC, double *ImC,  mwIndex *cindC, mwIndex *rindC)
{
  mwIndex iA, iB, jA, jB, jAmA, ck, Arow, Bend;
  mwIndex *Bptr;
  double Avalr, Avali;
  bool Acomp,Bcomp;
    
  if (ImA==NULL) {Acomp=false; Avali=0;} else Acomp=true;
  if (ImB==NULL) Bcomp=false; else Bcomp=true;

  ck=0;
  for (jA=0, jAmA=0; jA<nA; jA++, jAmA+=mA)
    for (jB=0, Bptr=cindB; jB<nB; jB++,Bptr++,*(++cindC)=ck)
      for (iA=0, Arow=0; iA<mA; iA++, Arow+=mB){
        Avalr=ReA[iA+jAmA]; if (Acomp) Avali=ImA[iA+jAmA];
        if (Avalr!=0 || Avali!=0)
          if (Bcomp)
            for (iB=*Bptr, Bend=*(Bptr+1); iB<Bend; iB++){
              *ReC++=Avalr*ReB[iB]-Avali*ImB[iB];
              *ImC++=Avali*ReB[iB]+Avalr*ImB[iB];
              *rindC++ = Arow+rindB[iB];
              ck++;
            }
          else
            for (iB=*Bptr, Bend=*(Bptr+1); iB<Bend; iB++){
              *ReC++=Avalr*ReB[iB];
              *ImC++=Avali*ReB[iB];
              *rindC++ = Arow+rindB[iB];
              ck++;
            }
      }
}


void kronssc(double *ReA, double *ImA, mwIndex *cindA, mwIndex *rindA, mwSize mA, mwSize nA,
             double *ReB, double *ImB, mwIndex *cindB, mwIndex *rindB, mwSize mB, mwSize nB, 
             double *ReC, double *ImC, mwIndex *cindC, mwIndex *rindC)
{
  mwIndex iA, iB, jA, jB, ck, Arow, Aend, Bend;
  mwIndex *Bptr;
  double Avalr, Avali;
  bool Acomp,Bcomp;
    
  if (ImA==NULL) {Acomp=false; Avali=0;} else Acomp=true;
  if (ImB==NULL) Bcomp=false; else Bcomp=true;

  ck=0;
  for (jA=0; jA<nA; jA++,cindA++)
    for (jB=0, Bptr=cindB; jB<nB; jB++,Bptr++,*(++cindC)=ck)
      for (iA=*cindA, Aend=*(cindA+1); iA<Aend; iA++){
        Avalr=ReA[iA]; if (Acomp) Avali=ImA[iA];
        Arow=rindA[iA]*mB;
        if (Bcomp)
          for (iB=*Bptr, Bend=*(Bptr+1); iB<Bend; iB++){
            *ReC++=Avalr*ReB[iB]-Avali*ImB[iB];
            *ImC++=Avali*ReB[iB]+Avalr*ImB[iB];
            *rindC++ = Arow+rindB[iB];
            ck++;
          }
        else
          for (iB=*Bptr, Bend=*(Bptr+1); iB<Bend; iB++){
            *ReC++=Avalr*ReB[iB];
            *ImC++=Avali*ReB[iB];
            *rindC++ = Arow+rindB[iB];
            ck++;
          }
      }
}



double *kronff(double *A, mwSize mA, mwSize nA, 
               double *B, mwSize mB, mwSize nB, 
               double *C)
{
  double *aptr, *bptr, *cptr;
  double Aval;
  mwIndex iA, iB, jA, jB;
  aptr=A;
  cptr=C;
  /* This is not needed because mxCreateDoubleMatrix does this */
  /*memset(C,0,mA*mB*nA*nB*sizeof(double));*/
  for (jA=0; jA<nA; jA++, aptr+=mA){
    for (jB=0, bptr=B; jB<nB; jB++, bptr+=mB){
      for (iA=0; iA<mA;){
        Aval=aptr[iA++];
        if (Aval!=0)
          for (iB=0; iB<mB;) *cptr++=Aval*bptr[iB++];
        else cptr+=mB;
      }
    }
  }
  return(C);
}

void kronsf(double *A, mwIndex *cindA, mwIndex *rindA, mwSize mA, mwSize nA, 
            double *B,                                 mwSize mB, mwSize nB, 
            double *C, mwIndex *cindC, mwIndex *rindC)
{
  mwIndex iA, iB, jA, jB, jBmB, ij, ck, Arow, Aend;
  double Aval;

  ck=0;
  for (jA=0; jA<nA; jA++,cindA++)
    for (jB=0, jBmB=0; jB<nB; jB++, jBmB+=mB, *(++cindC)=ck)
      for (iA=*cindA,Aend=*(cindA+1); iA<Aend; iA++){
        Aval=A[iA]; Arow=rindA[iA]*mB;
        for (iB=0, ij=jBmB; iB<mB; iB++, ij++)
          if (B[ij]!=0){
            *C++=Aval*B[ij];
            *rindC++ = Arow+iB;
            ck++;
          }
      }
}

void kronfs(double *A,                                 mwSize mA, mwSize nA,
            double *B, mwIndex *cindB, mwIndex *rindB, mwSize mB, mwSize nB,
            double *C, mwIndex *cindC, mwIndex *rindC)
{
  mwIndex iA, iB, jA, jB, jAmA, ck, Arow, Bend;
  mwIndex *Bptr;
  double Aval;

  ck=0;
  for (jA=0, jAmA=0; jA<nA; jA++, jAmA+=mA)
    for (jB=0, Bptr=cindB; jB<nB; jB++,Bptr++,*(++cindC)=ck)
      for (iA=0, Arow=0; iA<mA; iA++, Arow+=mB){
        Aval=A[iA+jAmA]; 
        if (Aval!=0)
          for (iB=*Bptr, Bend=*(Bptr+1); iB<Bend; iB++){
            *C++=Aval*B[iB];
            *rindC++ = Arow+rindB[iB];
            ck++;
          }   
      }
}


void kronss(double *A, mwIndex *cindA, mwIndex *rindA, mwSize mA, mwSize nA,
            double *B, mwIndex *cindB, mwIndex *rindB, mwSize mB, mwSize nB, 
            double *C, mwIndex *cindC, mwIndex *rindC)
{
  mwIndex iA, iB, jA, jB, ck, Arow, Aend, Bend;
  mwIndex *Bptr;
  double Aval;
  ck=0;
  for (jA=0; jA<nA; jA++,cindA++)
    for (jB=0, Bptr=cindB; jB<nB; jB++,Bptr++,*(++cindC)=ck)
      for (iA=*cindA, Aend=*(cindA+1); iA<Aend; iA++){
        Aval=A[iA]; Arow=rindA[iA]*mB;
        for (iB=*Bptr, Bend=*(Bptr+1); iB<Bend; iB++){
          *C++=Aval*B[iB];
          *rindC++ = Arow+rindB[iB];
          ck++;
        }
      }
}

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{  double *A, *B, *C;
   mwIndex *cindA, *cindB, *cindC, *rindA, *rindB, *rindC;
   mwSize mA, mB, nA, nB, nnzA, nnzB, mn;
   mwIndex i;

   if (nrhs!=2)
      mexErrMsgTxt("Two parameters must be passed");
   if (nlhs>1)
      mexErrMsgTxt("Only one output is created");
   if (!mxIsDouble(prhs[0]))
      mexErrMsgTxt("Inputs must have double data type"); 
   if (!mxIsDouble(prhs[1]))
      mexErrMsgTxt("Inputs must have double data type");
   
   A=mxGetPr(prhs[0]);
   mA=mxGetM(prhs[0]); 
   nA=mxGetN(prhs[0]);
   B=mxGetPr(prhs[1]);
   mB=mxGetM(prhs[1]);  
   nB=mxGetN(prhs[1]);   
   
   if (mA==0 || nA==0){
     plhs[0]=mxDuplicateArray(prhs[1]);
     return;
   } 
   if (mB==0 || nB==0){
     plhs[0]=mxDuplicateArray(prhs[0]);
     return;
   } 
   
   if (mxIsSparse(prhs[0]))
   {
     cindA=mxGetJc(prhs[0]);
     rindA=mxGetIr(prhs[0]);
     nnzA=cindA[nA];
     if (mxIsSparse(prhs[1]))
     {
       cindB=mxGetJc(prhs[1]);
       rindB=mxGetIr(prhs[1]);
       nnzB=cindB[nB];
       if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])){
         double *ImA, *ImB, *ImC;
           plhs[0]=mxCreateSparse(mA*mB,nA*nB,nnzA*nnzB,mxCOMPLEX);
           C=mxGetPr(plhs[0]);
           cindC=mxGetJc(plhs[0]);
           rindC=mxGetIr(plhs[0]);
           ImA=mxGetPi(prhs[0]);
           ImB=mxGetPi(prhs[1]);
           ImC=mxGetPi(plhs[0]);
           kronssc(A,ImA,cindA,rindA,mA,nA,B,ImB,cindB,rindB,mB,nB,C,ImC,cindC,rindC);
       }
       else{
         plhs[0]=mxCreateSparse(mA*mB,nA*nB,nnzA*nnzB,mxREAL);
         C=mxGetPr(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kronss(A,cindA,rindA,mA,nA,B,cindB,rindB,mB,nB,C,cindC,rindC);
       }
     }
     else
       if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])){
         double *ImA, *ImB, *ImC;
         ImA=mxGetPi(prhs[0]);
         ImB=mxGetPi(prhs[1]);
         nnzB=0;
         if (mxIsComplex(prhs[1])){
            for (i=0, mn=mB*nB; i<mn; i++){ if (B[i]!=0.0 || ImB[i]!=0.0) nnzB++;}
         }
         else{ 
            for (i=0, mn=mB*nB; i<mn; i++){if (B[i]!=0.0) nnzB++;}
         }
         plhs[0]=mxCreateSparse(mA*mB,nA*nB,nnzA*nnzB,mxCOMPLEX);
         C=mxGetPr(plhs[0]);
         ImC=mxGetPi(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kronsfc(A,ImA,cindA,rindA,mA,nA,B,ImB,mB,nB,C,ImC,cindC,rindC);
       }
       else{
         nnzB=0;
         for (i=0, mn=mB*nB; i<mn; i++) if (B[i]!=0.0) nnzB++;
         plhs[0]=mxCreateSparse(mA*mB,nA*nB,nnzA*nnzB,mxREAL);
         C=mxGetPr(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kronsf(A,cindA,rindA,mA,nA,B,mB,nB,C,cindC,rindC);     
       }
   }
   else{
     if (mxIsSparse(prhs[1])){
       cindB=mxGetJc(prhs[1]);
       rindB=mxGetIr(prhs[1]);
       nnzB=cindB[nB];
       nnzA=0;
       if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])){
         double *ImA, *ImB, *ImC;
         ImA=mxGetPi(prhs[0]);
         ImB=mxGetPi(prhs[1]);
         if (mxIsComplex(prhs[0])){
           for (i=0, mn=mA*nA; i<mn; i++){ if (A[i]!=0.0 || ImA[i]!=0.0) nnzA++;}
         }
         else{
           for (i=0, mn=mA*nA; i<mn; i++) if (A[i]!=0.0) nnzA++;
         }
         plhs[0]=mxCreateSparse(mA*mB,nA*nB,nnzA*nnzB,mxCOMPLEX);
         C=mxGetPr(plhs[0]);
         ImC=mxGetPi(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kronfsc(A,ImA,mA,nA,B,ImB,cindB,rindB,mB,nB,C,ImC,cindC,rindC);
       }
       else{
         for (i=0, mn=mA*nA; i<mn; i++) if (A[i]!=0.0) nnzA++;
         plhs[0]=mxCreateSparse(mA*mB,nA*nB,nnzA*nnzB,mxREAL);
         C=mxGetPr(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kronfs(A,mA,nA,B,cindB,rindB,mB,nB,C,cindC,rindC);
       }
     }
     else
     {
       if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])){
         double *ImA, *ImB, *ImC;
           plhs[0]=mxCreateDoubleMatrix(mA*mB,nA*nB,mxCOMPLEX);
           C=mxGetPr(plhs[0]);
           ImA=mxGetPi(prhs[0]);
           ImB=mxGetPi(prhs[1]);
           ImC=mxGetPi(plhs[0]);
           kronffc(A,ImA,mA,nA,B,ImB,mB,nB,C,ImC);
       }
       else{
         plhs[0]=mxCreateDoubleMatrix(mA*mB,nA*nB,mxREAL);
         C=mxGetPr(plhs[0]);
         kronff(A,mA,nA,B,mB,nB,C);
       }
     }
   }
}
