#include "mex.h"
#include <math.h>

#ifdef old32
#define mwSize int
#define mwIndex int
#endif

/* 
KRONCOL column Kronecker product
 C=kroncol(A,B);
 A and B are maxn and mbxn, C is ma*mbxn
*/


void kronffc(double *ReA, double *ImA, mwSize mA, 
             double *ReB, double *ImB, mwSize mB, 
             double *ReC, double *ImC, mwSize n)
{
  mwIndex iA, iB, j;
  double Avalr, Avali;
  bool Acomp,Bcomp;
    
  if (ImA==NULL) {Acomp=false; Avali=0;} else Acomp=true;
  if (ImB==NULL) Bcomp=false; else Bcomp=true;

  for (j=0; j<n; j++){
    for (iA=0; iA<mA; iA++){
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
    ReA+=mA;
    if (Acomp) ImA+=mA;
    ReB+=mB;
    if (Bcomp) ImB+=mB;
  }
}

void kronsfc(double *ReA, double *ImA,  mwIndex *cindA, mwIndex *rindA, mwSize mA, 
             double *ReB, double *ImB,                                  mwSize mB, 
             double *ReC, double *ImC,  mwIndex *cindC, mwIndex *rindC, mwSize n)
{
  mwIndex iA, iB, j, ck, Arow, Aend;
  double Avalr, Avali;
  bool Acomp,Bcomp;
    
  if (ImA==NULL) {Acomp=false; Avali=0;} else Acomp=true;
  if (ImB==NULL) Bcomp=false; else Bcomp=true;

  ck=0;
  for (j=0; j<n; j++,cindA++,*(++cindC)=ck){
      for (iA=*cindA,Aend=*(cindA+1); iA<Aend; iA++){
        Avalr=ReA[iA]; if (Acomp) Avali=ImA[iA]; 
        if (Avalr!=0 || Avali!=0){
          Arow=rindA[iA]*mB;
          if (Bcomp)
            for (iB=0; iB<mB; iB++){
              if (ReB[iB]!=0 || ImB[iB]!=0){
                *ReC++=Avalr*ReB[iB]-Avali*ImB[iB];
                *ImC++=Avalr*ImB[iB]+Avali*ReB[iB];
                *rindC++ = Arow+iB;
                ck++;
              }
            }
          else 
            for (iB=0; iB<mB; iB++){
              if (ReB[iB]!=0){
                *ReC++=Avalr*ReB[iB];
                *ImC++=Avali*ReB[iB];
                *rindC++ = Arow+iB;
                ck++;
              }
            }
        }
      }
    ReB+=mB;
    if (Bcomp) ImB+=mB;
  }
}

void kronfsc(double *ReA, double *ImA,                                  mwSize mA,
             double *ReB, double *ImB,  mwIndex *cindB, mwIndex *rindB, mwSize mB,
             double *ReC, double *ImC,  mwIndex *cindC, mwIndex *rindC, mwSize n)
{
  mwIndex iA, iB, j, ck, Arow, Bend;
  double Avalr, Avali;
  bool Acomp,Bcomp;
    
  if (ImA==NULL) {Acomp=false; Avali=0;} else Acomp=true;
  if (ImB==NULL) Bcomp=false; else Bcomp=true;

  ck=0;
  for (j=0; j<n; j++,*(++cindC)=ck,cindB++){
    for (iA=0, Arow=0; iA<mA; iA++, Arow+=mB){
      Avalr=ReA[iA]; if (Acomp) Avali=ImA[iA];
      if (Avalr!=0 || Avali!=0)
        if (Bcomp)
          for (iB=*cindB, Bend=*(cindB+1); iB<Bend; iB++){
            *ReC++=Avalr*ReB[iB]-Avali*ImB[iB];
            *ImC++=Avali*ReB[iB]+Avalr*ImB[iB];
            *rindC++ = Arow+rindB[iB];
            ck++;
          }
        else
          for (iB=*cindB, Bend=*(cindB+1); iB<Bend; iB++){
            *ReC++=Avalr*ReB[iB];
            *ImC++=Avali*ReB[iB];
            *rindC++ = Arow+rindB[iB];
            ck++;
          }
    }
    ReA+=mA;
    if (Acomp) ImA+=mA;
  }
}


void kronssc(double *ReA, double *ImA, mwIndex *cindA, mwIndex *rindA, mwSize mA,
             double *ReB, double *ImB, mwIndex *cindB, mwIndex *rindB, mwSize mB, 
             double *ReC, double *ImC, mwIndex *cindC, mwIndex *rindC, mwSize n)
{
  mwIndex iA, iB, j, ck, Arow, Aend, Bend;
  double Avalr, Avali;
  bool Acomp,Bcomp;
    
  if (ImA==NULL) {Acomp=false; Avali=0;} else Acomp=true;
  if (ImB==NULL) Bcomp=false; else Bcomp=true;

  ck=0;
  for (j=0; j<n; j++,cindA++,cindB++,*(++cindC)=ck)
    for (iA=*cindA, Aend=*(cindA+1); iA<Aend; iA++){
      Avalr=ReA[iA]; if (Acomp) Avali=ImA[iA];
      Arow=rindA[iA]*mB;
      if (Bcomp)
        for (iB=*cindB, Bend=*(cindB+1); iB<Bend; iB++){
          *ReC++=Avalr*ReB[iB]-Avali*ImB[iB];
          *ImC++=Avali*ReB[iB]+Avalr*ImB[iB];
          *rindC++ = Arow+rindB[iB];
          ck++;
        }
      else
        for (iB=*cindB, Bend=*(cindB+1); iB<Bend; iB++){
          *ReC++=Avalr*ReB[iB];
          *ImC++=Avali*ReB[iB];
          *rindC++ = Arow+rindB[iB];
          ck++;
        }
    }
}



void kroncolff(double *A, mwSize mA,  
               double *B, mwSize mB, 
               double *C, mwSize n)
{
  double *cptr,*Aend,*Bptr,*Bend;
  double Aval;
  mwIndex j;
  cptr=C;
  Aend=A;
  Bend=B;
  for (j=0; j<n; j++){
    B=Bend; Bend=B+mB;
    for (Aend+=mA; A<Aend;){
      Aval=*A++;
      if (Aval!=0)
        for (Bptr=B; Bptr<Bend;) *cptr++ = Aval * *Bptr++;
      else cptr+=mB;
    }
  }
}


void kroncolfs(double *A,                                 mwSize mA, 
               double *B, mwIndex *cindB, mwIndex *rindB, mwSize mB,
               double *C, mwIndex *cindC, mwIndex *rindC, mwSize n)
{
  mwIndex j, nj, Arow, Aend, *rptr;
  double Aval, *Bptr, *Bend, *Cptr;
  Aend=mA*mB;
  Bend=B;
  Cptr=C;
  for (j=0; j<n; j++, rindB+=nj, cindB++){
    nj=(*(cindB+1)-*(cindB));
    B=Bend; Bend+=nj; 
    for (Arow=0; Arow<Aend; Arow+=mB){
      Aval=*A++; 
      if (Aval!=0)
        for (Bptr=B, rptr=rindB; Bptr<Bend; ){
          *Cptr++  = Aval * *Bptr++;
          *rindC++ = Arow + *rptr++;
        }   
    }
    *(++cindC)=(mwIndex)(Cptr-C);
  }
}

void kroncolsf(double *A, mwIndex *cindA, mwIndex *rindA, mwSize mA,  
               double *B,                                 mwSize mB, 
               double *C, mwIndex *cindC, mwIndex *rindC, mwSize n)
{
  mwIndex j, Crow, rowend;
  double Aval, *Aptr, *Aend, *Bptr, *Cptr;
  Aptr=A;
  Cptr=C;
  for (j=0; j<n; j++, B+=mB){
    for (Aend=A+*(++cindA); Aptr<Aend; ){
      Aval=*Aptr++; Crow = mB * *rindA++;
      for (Bptr=B, rowend=Crow+mB; Crow<rowend;){
        if (*Bptr!=0){
          *Cptr++   = Aval * *Bptr++;
          *rindC++  = Crow++;
        }
        else{
          Bptr++;
          Crow++;
        }
      }
    }
    *(++cindC)=(mwIndex)(Cptr-C);
  }
}

void kroncolss(double *A, mwIndex *cindA, mwIndex *rindA, mwSize mA, 
               double *B, mwIndex *cindB, mwIndex *rindB, mwSize mB, 
               double *C, mwIndex *cindC, mwIndex *rindC, mwSize n)
{
  mwIndex j, nj, Arow, *rptr;
  double Aval, *Aptr, *Aend, *Bptr, *Bend, *Cptr;
  Aptr=A;
  Bend=B;
  Cptr=C;
  for (j=0; j<n; j++, cindB++, rindB+=nj){
    nj=(*(cindB+1)-*(cindB));
    B=Bend; Bend += nj;
    for (Aend = A+*(++cindA); Aptr<Aend; ){
      Aval = *Aptr++; Arow = mB * *rindA++; 
      for (Bptr=B, rptr=rindB; Bptr<Bend;){
        *Cptr++  = Aval * *Bptr++;
        *rindC++ = Arow + *rptr++;
      }
    }
    *(++cindC)=(mwIndex)(Cptr-C);
  }
}
      
      
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{  double *A, *B, *C;
   mwIndex *cindA, *cindB, *cindC, *rindA, *rindB, *rindC;
   mwSize mA, mB, n, nB, nnz;
   mwIndex i;

   if (nrhs!=2)
      mexErrMsgTxt("Two parameters must be passed");
   if (nlhs>1)
      mexErrMsgTxt("Only one output is created");
   if (!mxIsDouble(prhs[0]))
      mexErrMsgTxt("Inputs must be matrices"); 
   if (!mxIsDouble(prhs[1]))
      mexErrMsgTxt("Inputs must be matrices");
  
   A=mxGetPr(prhs[0]);
   mA=mxGetM(prhs[0]); 
   B=mxGetPr(prhs[1]);
   mB=mxGetM(prhs[1]);  
   n=mxGetN(prhs[0]);
   nB=mxGetN(prhs[1]);
   
   
   if (mA==0 || n==0){
     plhs[0]=mxDuplicateArray(prhs[1]);
     return;
   } 
   if (mB==0 || nB==0){
     plhs[0]=mxDuplicateArray(prhs[0]);
     return;
   } 
   
   if (n!=nB)
      mexErrMsgTxt("Matrices must have the same number of columns");
   
   
   if (mxIsSparse(prhs[0])){
     cindA=mxGetJc(prhs[0]);
     rindA=mxGetIr(prhs[0]);
     if (mxIsSparse(prhs[1])){
       cindB=mxGetJc(prhs[1]);
       rindB=mxGetIr(prhs[1]);
       nnz=0;
       for (i=1;i<=n;i++) nnz += (cindA[i]-cindA[i-1])*(cindB[i]-cindB[i-1]);
       if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])){
         double *ImA, *ImB, *ImC;
           plhs[0]=mxCreateSparse(mA*mB,n,nnz,mxCOMPLEX);
           C=mxGetPr(plhs[0]);
           cindC=mxGetJc(plhs[0]);
           rindC=mxGetIr(plhs[0]);
           ImA=mxGetPi(prhs[0]);
           ImB=mxGetPi(prhs[1]);
           ImC=mxGetPi(plhs[0]);
           kronssc(A,ImA,cindA,rindA,mA,B,ImB,cindB,rindB,mB,C,ImC,cindC,rindC,n);
       }
       else{
         plhs[0]=mxCreateSparse(mA*mB,n,nnz,mxREAL);
         C=mxGetPr(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kroncolss(A,cindA,rindA,mA,B,cindB,rindB,mB,C,cindC,rindC,n);
       }
     }
     else{
       nnz=cindA[n]*mB;
       if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])){
         double *ImA, *ImB, *ImC;
         ImA=mxGetPi(prhs[0]);
         ImB=mxGetPi(prhs[1]);
         plhs[0]=mxCreateSparse(mA*mB,n,nnz,mxCOMPLEX);
         C=mxGetPr(plhs[0]);
         ImC=mxGetPi(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kronsfc(A,ImA,cindA,rindA,mA,B,ImB,mB,C,ImC,cindC,rindC,n);
       }
       else{
         plhs[0]=mxCreateSparse(mA*mB,n,nnz,mxREAL);
         C=mxGetPr(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kroncolsf(A,cindA,rindA,mA,B,mB,C,cindC,rindC,n);    
       }
     }
   }
   else{
     if (mxIsSparse(prhs[1])){
       cindB=mxGetJc(prhs[1]);
       rindB=mxGetIr(prhs[1]);
       nnz=cindB[n]*mA;
       if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])){
         double *ImA, *ImB, *ImC;
         ImA=mxGetPi(prhs[0]);
         ImB=mxGetPi(prhs[1]);
         plhs[0]=mxCreateSparse(mA*mB,n,nnz,mxCOMPLEX);
         C=mxGetPr(plhs[0]);
         ImC=mxGetPi(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kronfsc(A,ImA,mA,B,ImB,cindB,rindB,mB,C,ImC,cindC,rindC,n);
       }
       else{
         plhs[0]=mxCreateSparse(mA*mB,n,nnz,mxREAL);
         C=mxGetPr(plhs[0]);
         cindC=mxGetJc(plhs[0]);
         rindC=mxGetIr(plhs[0]);
         kroncolfs(A,mA,B,cindB,rindB,mB,C,cindC,rindC,n);
       }
     }
     else{
       if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])){
         double *ImA, *ImB, *ImC;
           plhs[0]=mxCreateDoubleMatrix(mA*mB,n,mxCOMPLEX);
           C=mxGetPr(plhs[0]);
           ImA=mxGetPi(prhs[0]);
           ImB=mxGetPi(prhs[1]);
           ImC=mxGetPi(plhs[0]);
           kronffc(A,ImA,mA,B,ImB,mB,C,ImC,n);
       }
       else{
         plhs[0]=mxCreateDoubleMatrix(mA*mB,n,mxREAL);
         C=mxGetPr(plhs[0]);
         kroncolff(A,mA,B,mB,C,n);
       }
     }
   }
}
