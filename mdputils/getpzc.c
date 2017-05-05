#include "mex.h"
#include <math.h>

// include tmwtypes.h and use int16_T, uint16_T, int32_T, uint32_T, real_T etc. 


//
// This is the basic MEX file used by catcountP to
// compute probabilities. It's calling syntax is
//    Pz=getpzc(S,xind,nxind,tab,factor,logp,SXj);
// where
//    S      : a grid of state values
//    xind   : a cell array of index vectors to extract columns of S
//    nxind  : a vector with the number of elements in each of the xind index vectors
//    tab    : a table of multiset coefficients for determining state locations
//    factor : vector of precomputed factors used to compute multinomial probabilities
//    logp   : nxm matrix of log probability values
//    SXj    : m-vector of state/action combinations
//
// This procedure has no checks and should generally not be called directly. 
// It is called by catcountP
// 

double *logp, *factor, ninf;

unsigned int *tab, *tabstart, indstart, p2;

//  quickSort
//
//  This public-domain C implementation by Darel Rex Finley.
//
//  * Returns YES if sort was successful, or NO if the nested
//    pivots went too deep, in which case your array will have
//    been re-ordered, but probably not sorted correctly.
//
//  * This function assumes it is called with valid parameters.
//
//  * Example calls:
//    quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7

// http://www.alienryderflex.com/quicksort/

// modified to also return the index for the sort order,
// sort in descending order and input unsigned integers to be sorted

void quickSort(unsigned int *arr, unsigned int *ind, int elements) {
  #define  MAX_LEVELS  300
  unsigned int  piv, pivind;
  int i=0, beg[MAX_LEVELS], end[MAX_LEVELS], L, R, swap;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=arr[L]; pivind=ind[L];
      while (L<R) {
        while (arr[R]<=piv && L<R) R--; if (L<R) {ind[L]=ind[R];arr[L++]=arr[R];}
        while (arr[L]>=piv && L<R) L++; if (L<R) {ind[R]=ind[L];arr[R--]=arr[L];}
      }
      arr[L]=piv;  ind[L]=pivind;
      beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
      if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; 
      }
    }
    else i--;
  }
}


void insertionSort(unsigned int *a, unsigned int *ind, int n){
  int i, j;
  unsigned int ax, ix;

  for (i=1; i < n; i++) {
    ax = a[i]; ix = ind[i];
    j = i;
    while ((j >= 0) && (a[j-1] < ax)) {
      a[j]   = a[j-1];
      ind[j] = ind[j-1];
      j--;
    }
    a[j]   = ax;
    ind[j] = ix;
  }
}

unsigned int getind(unsigned int *x, unsigned int *y){
unsigned int  ind, *tabptr;
  tabptr = tabstart;
  ind    = indstart;
  while (tabptr>=tab){
    tabptr -= (unsigned int) (*x++ + *y++);  
    ind    -= *tabptr;
    tabptr -= p2;
  }
  return(ind);
}


double multinomialOLD(
  unsigned int *x, unsigned int n, unsigned int Sji, unsigned int Sii)
{
  unsigned int xi, xend;
  double p, *lp;
  mwIndex i;
  
  lp=logp+Sii*(n+1);
  p=factor[Sji];
  xend=Sji;
  for (i=0;i<n;i++){
    xi = *x++;
    xend -= xi;
    if (xi==0) p -= factor[0];               // avoid 0 times -inf problems
    else       p += lp[i]*xi - factor[xi];
  }
  if (xend==0) p -= factor[0];
  else         p += lp[n]*xend - factor[xend];
  return(exp(p));
}

double multinomial(
  unsigned int *x, unsigned int n, unsigned int Sji, unsigned int Sii)
{
  unsigned int xi, xend;
  double p, *lp;
  mwIndex i;
  
  lp=logp+Sii*(n+1);
  p=factor[Sji];
  xend=Sji;
  for (i=0;i<n;i++){
    xi = *x++;
    xend -= xi;
    if (xi==0) p -= factor[0];               // avoid 0 times -inf problems
    else{
      if (lp[i]==ninf) return(0);
      p += lp[i]*xi - factor[xi];
    }
  }
  if (xend==0) p -= factor[0];
  else         p += lp[n]*xend - factor[xend];
  return(exp(p));
}

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double *p, *pz, pxj, *ptemp;
  mwSize n, n1, n2, m, nm1, pn;
  mwIndex i, j, ii;
  unsigned int *z, *x, *S, B1, Sji, Sii, indpij, lastind, *nxind, *zind, *xind;
  unsigned int *Sj, *Si, mm, nnz;
  bool *skip;

  // Error checking on inputs  
  if (nrhs<7 || nrhs>7) mexErrMsgTxt("Incorrect number of input arguments.");
  
  
  //Pj=getpzc2(S,xind,nxind,tab,factor,logp,uint32(Xj))
  
  n=mxGetM(prhs[0]);
  nm1=n-1;
  m=mxGetM(prhs[3]); 
  S      = mxGetData(prhs[0]);
  tab    = mxGetData(prhs[3]);
  factor = mxGetPr(prhs[4]);
  logp   = mxGetPr(prhs[5]);
  
  //need to put in check that X is uint32 and, if not, convert it
  Sj     = mxGetData(prhs[6]);
  
  mm=mxGetM(prhs[6]);
  p2=mxGetM(prhs[3]);
  nxind=mxGetData(prhs[2]);
  pn=nxind[mxGetNumberOfElements(prhs[2])-1];
  
  Si=mxCalloc(mm,sizeof(unsigned int));
  pz=mxCalloc(pn,sizeof(double));
  skip=mxCalloc(pn,sizeof(bool));
  
  for (i=0; i<mm; i++) Si[i]=i;
  // Use insertion sort when mm is small, otherwise use quick sort
  if (mm>16)
    quickSort(Sj, Si, mm); 
  else
    insertionSort(Sj, Si, mm);
  plhs[0]=mxCreateDoubleMatrix(pn,1,mxREAL);  
  p=mxGetPr(plhs[0]);
  ninf=-mxGetInf();
                
  Sji=Sj[0];
  Sii=Si[0];
  xind=mxGetData(mxGetCell(prhs[1],Sji));
  n1=nxind[Sji];
  nnz=0; for (i=0; i<mm; i++) if (Sj[i]>0) nnz++;  else break;
  // if nnz is an even # of non-zeros swap p and pz
  // this avoid copying p into pz and ensures that the last time though
  // the loop p is the pointer to the output vector
  if ((nnz/2)*2==nnz){ptemp=p; p=pz; pz=ptemp;}  
  for (i=0; i<n1; i++) p[i] = multinomial(S+xind[i]*n, n, Sji, Sii);    
  B1=Sji;          
  for (ii=1;ii<nnz;ii++){   // only process non-zero values of Xj
    Sji=Sj[ii];
    Sii=Si[ii];
    n1=nxind[B1];
    n2=nxind[Sji];
    zind=mxGetData(mxGetCell(prhs[1],B1));
    xind=mxGetData(mxGetCell(prhs[1],Sji));
    ptemp=p; p=pz; pz=ptemp;                   // swap p pointers
    ptemp=p+n1; while (p<ptemp) *(--ptemp)=0;  // set p to 0
    B1 += Sji;
    tabstart=tab+p2*(n-1)+B1;
    indstart=*(tabstart+1) - 1;
    // process the first element
    x=S+xind[0]*n;
    pxj = multinomial(x, n, Sji, Sii);
    p[0] = pz[0]*pxj;
    indpij=0;
    skip[0]=false;
    // process the next n1-1 elements to determine the skip pattern
    for (i=1; i<n1; i++){
      lastind=indpij;
      z = S+zind[i]*n;
      indpij=getind(x,z);
      p[indpij] = pz[i]*pxj;
      if (indpij==lastind+1) skip[i]=true;
      else                   skip[i]=false;
    }
    // loop over the reamining n2 x values
    for (j=1; j<n2; j++){
      x = S+xind[j]*n;
      pxj = multinomial(x, n, Sji, Sii);
      for (i=0; i<n1; i++){
        if (skip[i]) indpij++;
        else         indpij=getind(x,S+zind[i]*n);
        p[indpij] += pz[i]*pxj;
      }
    }
  }
  mxFree(skip);
  mxFree(pz);
  mxFree(Si);
}


