#include "mex.h"

/*
 * VestBMS_finalqtrapz.c
 *
 * multiplies three matrices elementwise and executes trapezoidal 
 * integration along second and third dimensions.
 *
 * This is a MEX-file for MATLAB.
 */


void qtrapzc(double *xpdf_a, double *xpdf_b, double *w, double *z, mwSize N, mwSize K)
{
    mwSize i,j,k,m;
    double sum,jsum,*w0;
    
    w0 = w;     /* copy base pointer */
    
    for (i=0; i<N; i++) {
        sum = 0.;
        w = w0;     /* reset pointer before inner loops */
        for (k=0; k<K; k++) {
            
            /* inner qtrapz loop */
            jsum = 0.5 * xpdf_a[i] * *(w++);
            for (j=1; j<K-1; j++) {
                jsum += xpdf_a[i+j*N] * *(w++);
            }
            jsum += 0.5 * xpdf_a[i+(K-1)*N] * *(w++);
            
            if (k == 0 || k == K-1) {
                sum += 0.5 * jsum * xpdf_b[i+k*N];
            }
            else {
                sum += jsum * xpdf_b[i+k*N];                
            }
        }
        *(z++) = sum;
    }
        
    /* printf("%d %d\n",(int)idx[0],(int)idx[1]); */
    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  const mwSize *dims_a,*dims_b,*dims_w;
  double *xpdf_a,*xpdf_b,*w,*z;
  int N,K;
  size_t n1,n2,n3;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs<2 || nrhs>3) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "Two or three inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "One output required.");
  
  /* check to make sure the first input argument is a scalar */
  /* if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:xNotScalar",
            "Input x must be a scalar.");
  } */
    
  /* Create a pointer to the first matrix (xpdf_vis, N-by-K-by-1) */
  xpdf_a = mxGetPr(prhs[0]);
  dims_a = mxGetDimensions(prhs[0]);    /*  get the dimensions */

  /* Create a pointer to the second matrix (xpdf_vest, N-by-1-by-K) */
  xpdf_b = mxGetPr(prhs[1]);
  dims_b = mxGetDimensions(prhs[1]);    /*  get the dimensions */
  
  /* Create a pointer to the third matrix (w, 1-by-K-by-K) */
  w = mxGetPr(prhs[2]);
  dims_w = mxGetDimensions(prhs[2]);    /*  get the dimensions */
    
  /* Check matrix sizes */
  N = dims_a[0];
  K = dims_a[1];
  if ( dims_b[0] != N || dims_b[1] != 1 || dims_b[2] != K ) {
      mexErrMsgIdAndTxt( "MATLAB:VestBMS_finalqtrapzc:dimensionMismatch",
            "Second matrix dimensions do not match the first input.");      
  }
  if ( dims_w[0] != 1 || dims_w[1] != K || dims_w[2] != K ) {
      mexErrMsgIdAndTxt( "MATLAB:VestBMS_finalqtrapzc:dimensionMismatch",
            "Third matrix dimensions do not match the first input.");      
  }

  /*  set the output pointer to the output matrix (N-by-1) */
  plhs[0] = mxCreateDoubleMatrix((mwSize) N, 1, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  qtrapzc(xpdf_a, xpdf_b, w, z, (mwSize) N,(mwSize) K);
  
}
