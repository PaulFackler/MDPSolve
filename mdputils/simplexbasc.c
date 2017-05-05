#include "mex.h"
#include <math.h>

/* Form an interpolation basis for a simplex */


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
// and sort in descending order and to be stable (preserve original order for ties)
void quicksort(double *arr, mwSignedIndex *ind, mwSize elements) {
  #define  MAX_LEVELS  300
  double piv;
  mwSignedIndex pivind;
  mwSignedIndex i=0, beg[MAX_LEVELS], end[MAX_LEVELS], L, R, swap;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=arr[L]; pivind=ind[L];
      while (L<R) {
        while ((arr[R]<piv || (arr[R]==piv && ind[R]>pivind)) && L<R) R--; 
        if (L<R) {ind[L]=ind[R];arr[L++]=arr[R];}
        while ((arr[L]>piv || (arr[L]==piv && ind[L]<pivind)) && L<R) L++; 
        if (L<R) {ind[R]=ind[L];arr[R--]=arr[L];}
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


// modified to also return the index for the sort order,
// and sort in descending order 
void quicksort2(mwSignedIndex *arr, double *ind, mwSize elements) {
  #define  MAX_LEVELS  300
  double pivind;
  mwSignedIndex piv;
  mwSignedIndex i=0, beg[MAX_LEVELS], end[MAX_LEVELS], L, R, swap;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=arr[L]; pivind=ind[L];
      while (L<R) {
        while (arr[R]>=piv && L<R) R--; if (L<R) {ind[L]=ind[R];arr[L++]=arr[R];}
        while (arr[L]<=piv && L<R) L++; if (L<R) {ind[R]=ind[L];arr[R--]=arr[L];}
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

// sorts a in descending order
// also returns the index of the sorted values
// ind should start as 1:n
void insertionsort(double *a, mwIndex *ind, mwIndex n){
  mwIndex i, j, j1;
  mwIndex ix;
  double ax;
  for (i=1; i < n; i++) {
    ax = a[i]; ix = ind[i];
    j = i; j1=j-1;
    while (a[j1] < ax) {
      a[j]   = a[j1];
      ind[j] = ind[j1];
      j=j1; 
      if (j1==0) break;
      else j1--;
    }
    a[j]   = ax;
    ind[j] = ix;
  }
}

// sorts a so its indices are in acsending order
// also returns the indices 
// note that a and ind are reverses here
void insertionsort2(mwIndex *a, double *ind, mwIndex n){
  mwIndex i, j, j1;
  mwIndex ax;
  double ix;
  for (i=1; i < n; i++) {
    ax = a[i]; ix = ind[i];
    j = i; j1 = j-1;
    while (a[j1] > ax) {
      a[j]   = a[j1];
      ind[j] = ind[j1];
      j=j1; 
      if (j1==0) break;
      else j1--;
    }
    a[j]   = ax;
    ind[j] = ix;
  }
}


 mwIndex *tab;
 mwSize q, q1, p, p1;
 
 // table of multiset coefficients
// table in column reversed order
mwIndex maketab( mwSize p1,  mwSize q1){
mwIndex i, j, *tabptr;
  if (q1==1) return(p1);
  tabptr=tab+(q1-1)*p1-1;
  tabptr[1]=1;
  for (i=2; i<=p1; i++) tabptr[i] = tabptr[i-1]+1;
  for (j=1; j<q1; j++){
    tabptr -= p1;
    tabptr[1] = 1;
    for (i=2; i<=p1; i++) tabptr[i] = tabptr[i-1] + tabptr[i+p1];
  }
  return(tabptr[p1]);
}

// find index of vertex v
 mwIndex getind(mwIndex *v){
 mwIndex  ind, *tabptr, i, iend;
 mwSignedIndex  eta;
  tabptr = tab;
  ind    = tab[p1];
  eta    = p;
  iend   = q1-1;
  for (i=0; i<iend; i++){
    eta    -= v[i] - v[i+1];
    if (eta<=0) return(ind);
    ind    -= tabptr[eta];
    tabptr += p1;
  }
  eta -= v[iend];
  if (eta>0) ind -= tabptr[eta];
  return(ind);
}


void simplexbas(double *xj, double *w, mwIndex *ir, mwIndex *v){
mwIndex i;
double wi;
  if (q1==1){    // faster operation when q=2 
    wi=xj[0];
    if (wi<0) *ir=0;
    else
      if (wi>=p) *ir=p-1;
      else *ir = floor(wi);
    wi -= *ir;
    (*ir)++; 
    ir[1] = (*ir)+1; 
    w[1]=wi;
    w[0]=1-wi;
  }
  else{
    i=q1;
    wi=0;
    ir++;
    while (i>0){
      i--; 
      wi += xj[i];   // reverse cumulative sum of xj
      if (wi<0) v[i]=0;
      else{ 
        if (wi>=p) v[i]=p-1;
        else       v[i]=floor(wi);
      }
      w[i] = wi - v[i];
      ir[i]=i;
    }
    //quicksort(w,ir,q1);
    insertionsort(w,ir,q1);    
    w[q1]=w[q1-1];
    for (i=q1-1; i>0; i--){w[i] = w[i-1]-w[i];}
    w[0]=1-w[0];
    if (w[0]<0 && w[0]>-1e-12) w[0]=0; // check for small negative numbers
    ir--;
    ir[0]=getind(v); 
    for (i=1; i<=q1;i++){
      v[ir[i]] ++;
      ir[i]=getind(v); 
    }
    //quicksort2(ir,w,q);  // sort by rows of each column so data is ready for sparse form
    insertionsort2(ir,w,q);  // sort by rows of each column so data is ready for sparse form
  }
}




void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  int ii;
  double *x, *xptr, *xj, C, *w, *wj;
  mwIndex *ir, *jc, *irj, i, j, *v;
  mwSize N, n;

  /* Error checking on inputs */  
  if (nrhs<3) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>4) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>2) mexErrMsgTxt("Too many output arguments.");
  for (ii=0; ii<nrhs; ii++){
    if (!mxIsDouble(prhs[ii]) || mxIsSparse(prhs[ii]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[ii]))
      mexErrMsgTxt("X must be real.");
  }
  if (mxGetNumberOfElements(prhs[1])>1)
      mexErrMsgTxt("Second input must be scalar");
  if (mxGetNumberOfElements(prhs[2])>1)
      mexErrMsgTxt("Third input must be scalar");
  if (nrhs>=4)
    if (mxGetNumberOfElements(prhs[3])>1)
       mexErrMsgTxt("Fourth input must be scalar");
  
  x=mxGetPr(prhs[0]);
  q=(mwSize) *mxGetPr(prhs[1]);
  if (q<2)
     mexErrMsgTxt("q must be at least 2");
  p=(mwSize) *mxGetPr(prhs[2]);
  if (nrhs>=4) C=*mxGetPr(prhs[3]);
  else         C=(double) p;
  
  q1=q-1;
  p1=p+1;
  N=mxGetM(prhs[0]);                    // x should be N x q (or N x q-1)
  if (mxGetN(prhs[0])>q || mxGetN(prhs[0])<q1)
       mexErrMsgTxt("Input data should have q columns");
  
  // create a table of multiset coefficients
  if (q1>1){
    tab=mxCalloc(p1*q1,sizeof(mwIndex));
    n=maketab(p1,q1);
    tab--;
  }
  // allocation workspace memory  
  xj=mxCalloc(q-1,sizeof(double));
  v=mxCalloc(q,sizeof(mwIndex));
  
  // allocate memory for outputs
  if (nlhs<2){
    if (q1==1) n=p1;
    plhs[0]=mxCreateSparse(n, N, q*N, mxREAL);
    w = mxGetPr(plhs[0]);
    ir = mxGetIr(plhs[0]);
    jc = mxGetJc(plhs[0]);
    for (j=0;j<=N;j++) jc[j]=j*q;
  }
  else{
    plhs[0]=mxCreateDoubleMatrix(q,N,mxREAL);
    w = mxGetPr(plhs[0]);
    if (sizeof(mwIndex)==4)
      plhs[1]=mxCreateNumericMatrix(q,N,mxUINT32_CLASS,mxREAL);
    else
      plhs[1]=mxCreateNumericMatrix(q,N,mxUINT64_CLASS,mxREAL);
    ir = mxGetData(plhs[1]);
  }
  
  // loop over the input data
  // If C==p there is no need to adjust the data
  wj=w;
  irj=ir;
  if (C==p){
    for (j=0; j<N; j++,wj+=q,irj+=q){
      xptr=x+j;
      for (i=0;i<q1;i++,xptr+=N) xj[i]=*xptr;
      simplexbas(xj,wj,irj,v);
    }
  }
  else{
    C=p/C;
    for (j=0; j<N; j++,wj+=q,irj+=q){
      xptr=x+j;
      for (i=0;i<q1;i++, xptr+=N) xj[i]=*xptr*C;
      simplexbas(xj,wj,irj,v);
    }
  }
  // convert to 1 based indexing
  if (nlhs<2){
    n=q*N; 
    for (i=0; i<n; i++){ 
      ir[i]--;
    }
  }
  mxFree(xj);
  mxFree(v);
  if (q1>1) { tab++; mxFree(tab);}
}

