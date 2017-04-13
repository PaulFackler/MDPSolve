#include "mex.h"
#include <math.h>

#ifdef old32
#define mwSize int
#define mwIndex int
#endif

/* 
C code to compute direct products (row-wise tensor products)
For real full or sparse matrices only
*/

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{  double *a, *b, *c;
   mwIndex i, j, k;
   mwSize m, an, bn;
   if (nrhs!=2){mexErrMsgTxt("Two parameters must be passed");}
   if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
     mexErrMsgTxt("Inputs must be full or sparse double matrices");

   if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
         mexErrMsgTxt("Complex matrices not supported in dprod");


   if (mxIsEmpty(prhs[0]))
      plhs[0]=mxDuplicateArray(prhs[1]); 
   else if (mxIsEmpty(prhs[1]))
      plhs[0]=mxDuplicateArray(prhs[0]);
   else
   {
     m=mxGetM(prhs[0]);
     if (mxGetM(prhs[1])!=m)
       {mexErrMsgTxt("Matrices must have the same number of rows");}
     /* SPARSE/SPARSE -> SPARSE */
     if (mxIsSparse(prhs[0]) && mxIsSparse(prhs[1]))
     {
       mwIndex *ai, *aj, *bi, *bj, *ci, *cj;
       mwIndex ak, anj, bk, bjend, ck;
       mwSize *acounts, *bcounts;
       a=mxGetPr(prhs[0]);
       ai=mxGetIr(prhs[0]);
       aj=mxGetJc(prhs[0]);
       an=mxGetN(prhs[0]);
       b=mxGetPr(prhs[1]);
       bi=mxGetIr(prhs[1]);
       bj=mxGetJc(prhs[1]);
       bn=mxGetN(prhs[1]);
       /* determine # of non-zeros */
       acounts=mxCalloc(2*m,sizeof(mwSize));
       bcounts=acounts+m;
       for (k=*(aj+an), i=0; i<k; i++) acounts[ai[i]]++; 
       for (k=*(bj+bn), i=0; i<k; i++) bcounts[bi[i]]++;
       for (k=0,i=0; i<m; i++) k+=acounts[i]*bcounts[i];
       mxFree(acounts);
       /* create output array */
       plhs[0]=mxCreateSparse(m,an*bn,k,mxREAL);
       ci=mxGetIr(plhs[0]);
       cj=mxGetJc(plhs[0]);
       c=mxGetPr(plhs[0]);
       /* main loop order: columns of a, columns of b, rows */
       *cj=0;
       ck=0;
       for (i=0; i<an; i++,a+=anj,ai+=anj){
         anj=aj[i+1]-aj[i];
         bjend=0;
         for (j=1; j<=bn; j++, *(++cj)=ck){
           bk=bjend; bjend=bj[j];
           for (ak=0; ak<anj && bk<bjend;){
             if (ai[ak]==bi[bk]){
               *c++=a[ak]*b[bk++];
               *ci++=ai[ak++];
               ck++;
             }
             else if (ai[ak]<bi[bk]) ak++;
             else bk++;
           } 
         }
       }
     } 
     /* DOUBLE/SPARSE -> SPARSE */
     else if (mxIsDouble(prhs[0]) && mxIsSparse(prhs[1]))
     {
       mwIndex *bi, *bj, *ci, *cj;
       mwSize ak, bk, ck, bjj;
       mwSize *acounts, *bcounts;
       a=mxGetPr(prhs[0]);
       an=mxGetN(prhs[0]);
       b=mxGetPr(prhs[1]);
       bi=mxGetIr(prhs[1]);
       bj=mxGetJc(prhs[1]);
       bn=mxGetN(prhs[1]);
       /* determine # of non-zeros */
       acounts=mxCalloc(2*m,sizeof(mwSize));
       bcounts=acounts+m;
       for (j=0; j<an; j++) for (i=0; i<m; i++) if (*a++!=0) acounts[i]++; 
       a=a-m*an;
       for (k=*(bj+bn), i=0; i<k; i++) bcounts[bi[i]]++;
       for (k=0,i=0; i<m; i++) k+=acounts[i]*bcounts[i];
       mxFree(acounts);

       /* create output array */
       plhs[0]=mxCreateSparse(m,an*bn,k,mxREAL);
       ci=mxGetIr(plhs[0]);
       cj=mxGetJc(plhs[0]);
       c=mxGetPr(plhs[0]);
       /* main loop order: columns of a, columns of b, rows */
       *cj++=0;
       ck=0;
       for (i=0; i<an; i++,a+=m){
         bk=0;
         for (j=1; j<=bn; j++){
           bjj=bj[j];
           for (;bk<bjj;bk++){
             ak=bi[bk];
             if (a[ak]!=0){
               *c++  = a[ak]*b[bk];
               *ci++ = ak;
               ck++;
             }
           } 
           *cj++=ck;
         }
       }
     }
     /* SPARSE/DOUBLE -> SPARSE */
     else if (mxIsSparse(prhs[0]) && mxIsDouble(prhs[1]))
     {
       mwIndex *ai, *aj, *ci, *cj;
       double *bj;
       mwIndex ak, bk, ck, akk;
       mwSize *acounts, *bcounts;
       a=mxGetPr(prhs[0]);
       ai=mxGetIr(prhs[0]);
       aj=mxGetJc(prhs[0]);
       an=mxGetN(prhs[0]);
       b=mxGetPr(prhs[1]);
       bn=mxGetN(prhs[1]);
       /* determine # of non-zeros */
       acounts=mxCalloc(2*m,sizeof(mwSize));
       bcounts=acounts+m;
       for (k=*(aj+an), i=0; i<k; i++) acounts[ai[i]]++;
       for (j=0; j<bn; j++) for (i=0; i<m; i++) if (*b++!=0) bcounts[i]++; 
       b=b-m*bn;
       for (k=0,i=0; i<m; i++) k+=acounts[i]*bcounts[i];
       mxFree(acounts);
       /* create output array */
       plhs[0]=mxCreateSparse(m,an*bn,k,mxREAL);
       ci=mxGetIr(plhs[0]);
       cj=mxGetJc(plhs[0]);
       c=mxGetPr(plhs[0]);
       /* main loop order: columns of a, columns of b, rows */
       *cj++=0;
       ck=0;
       for (i=0; i<an; i++, aj++){
         bj=b;
         for (j=1; j<=bn; j++, bj+=m){
           for (ak=*aj,akk=*(aj+1);ak<akk;ak++){
             bk=ai[ak];
             if (bj[bk]!=0){
               *c++  = a[ak]*bj[bk];
               *ci++ = bk;
               ck++;
             }
           } 
           *cj++=ck;
         }
       }
     } 
      /* DOUBLE/DOUBLE -> DOUBLE */
     else if (mxIsDouble(prhs[0]) && mxIsDouble(prhs[1])){
       mwSize mbn;
       double *aend;
       an=mxGetN(prhs[0]);
       bn=mxGetN(prhs[1]);
       mbn=bn*m;
       a=mxGetPr(prhs[0]);
       b=mxGetPr(prhs[1]);
       plhs[0]=mxCreateDoubleMatrix(m,an*bn,mxREAL);
       c=mxGetPr(plhs[0]);
       for (i = 0; i < an; i++, b-=mbn, a=aend)
         for (aend=a+m, j=0; j<bn; j++, a-=m)
           while (a<aend) *c++ = *a++ * *b++;
     }
     else mexErrMsgTxt("Unsupported data type passed to dprod");
   }
}
