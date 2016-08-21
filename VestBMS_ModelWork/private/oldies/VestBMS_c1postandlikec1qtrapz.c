#include "mex.h"
#include "float.h"

/*
 * VestBMS_c1postandlikec1qtrapz.c
 *
 * multiplies two matrices elementwise and compute posterior probability
 * of responding 'right' and returns likelihood.
 *
 * This is a MEX-file for MATLAB.
 */

void qtrapzc(double *x_a, double *x_b, double *z, double *z2, mwSize N, mwSize K)
{
    mwSize i,j,k,i0deg;
    double lsum,rsum,pmin,*x0_a,*x0_b;
    
    /* Index of 0 deg */
    i0deg = (N+1)/2; 
    
    /* Add some small probability to prevent 0 / 0 */
    pmin = 0.5 * DBL_MIN * (double) N;
    
    /* copy base pointers */    
    x0_a = x_a;
    x0_b = x_b;
    
    for (k=0; k<K; k++) {
        for (j=0; j<K; j++) {
            x_a = (x0_a + j*N);
            x_b = (x0_b + k*N);
            
            lsum = 0.5 * *(x_a++) * *(x_b++);
            for (i=1; i<i0deg-1; i++) {
                lsum += *(x_a++) * *(x_b++);
            }
            lsum += 0.5 * *(x_a) * *(x_b); /* last row, do not advance */
            
            rsum = 0.5 * *(x_a++) * *(x_b++);
            for (i=i0deg; i<N-1; i++) {
                rsum += *(x_a++) * *(x_b++);
            }
            rsum += 0.5 * *(x_a) * *(x_b);

            z[j*K+k] = (rsum + pmin) / (rsum + lsum + pmin);
            z2[j*K+k] = rsum + lsum;
        }
    }
        
    /* printf("%d %d\n",(int)idx[0],(int)idx[1]); */
    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  const mwSize *dims_a,*dims_b;
  double *x_a,*x_b,*w,*z,*z2;
  int N,K;
  size_t n1,n2,n3;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs<2 || nrhs>2) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "Two or three inputs required.");
  if(nlhs!=2) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "Two outputs required.");
  
  /* check to make sure the first input argument is a scalar */
  /* if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:xNotScalar",
            "Input x must be a scalar.");
  } */
    
  /* Create a pointer to the first matrix (postpdf_c2, N-by-1-by-K) */
  x_a = mxGetPr(prhs[0]);
  dims_a = mxGetDimensions(prhs[0]);    /*  get the dimensions */

  /* Create a pointer to the second matrix (like_vis, N-by-K) */
  x_b = mxGetPr(prhs[1]);
  dims_b = mxGetDimensions(prhs[1]);    /*  get the dimensions */
      
  /* Check matrix sizes */
  N = dims_a[0];
  K = dims_a[2];
  if ( dims_b[0] != N || dims_b[1] != K ) {
      mexErrMsgIdAndTxt( "MATLAB:VestBMS_likec1qtrapzc:dimensionMismatch",
            "Second matrix dimensions do not match the first input.");      
  }  
  /* printf("%d %d %d\n",dims_a[0],dims_a[1],dims_a[2]); */

  /*  set the output pointer to the first output matrix (K-by-K) */
  plhs[0] = mxCreateDoubleMatrix((mwSize) K, (mwSize) K, mxREAL);

  /*  set the output pointer to the second output matrix (K-by-K) */
  plhs[1] = mxCreateDoubleMatrix((mwSize) K, (mwSize) K, mxREAL);
  
  /* otherwise use mxCreateNumericMatrix(mwSize m, mwSize n, mxClassID classid, mxComplexity ComplexFlag); */
  
  /*  create a C pointer to a copy of the output matrices */
  z = mxGetPr(plhs[0]);
  z2 = mxGetPr(plhs[1]);
  
  /*  call the C subroutine */
  qtrapzc(x_a, x_b, z, z2, (mwSize) N,(mwSize) K);
  
}
