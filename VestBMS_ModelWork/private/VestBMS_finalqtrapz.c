#include "mex.h"
#include "math.h"
#include "matrix.h"

/*
 * VestBMS_finalqtrapz_mat_mex.c
 *
 *VESTBMS_FINALQTRAPZ Marginalize response probability over x_vis and x_vest
 *
 * ================ INPUT VARIABLES ====================
 * XPDF_VIS: p(x_vis|s). [S-by-K] (double)
 * XPDF_VEST: p(x_vest|s). [S-by-1-by-K] (double)
 * R: p(response|x_vis,x_vest). [1-by-K-by-K] (double)
 *
 * ================ OUTPUT VARIABLES ==================
 * PRMAT: p(r|s). [S-by-1] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 21-Aug-2016 with MEXXER v0.1 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 0

void VestBMS_finalqtrapz( double *prmat, double *xpdf_vis, double *xpdf_vest, double *R, mwSize K, mwSize S )
{
	
    mwSize i,j,k,m;
    double sum,jsum,*R0;
    
    R0 = R;     /* copy base pointer */
    
    for (i=0; i<S; i++) {
        sum = 0.;
        R = R0;     /* reset pointer before inner loops */
        for (k=0; k<K; k++) {
            
            /* inner qtrapz loop */
            jsum = 0.5 * xpdf_vis[i] * *(R++);
            for (j=1; j<K-1; j++) {
                jsum += xpdf_vis[i+j*S] * *(R++);
            }
            jsum += 0.5 * xpdf_vis[i+(K-1)*S] * *(R++);
            
            if (k == 0 || k == K-1) {
                sum += 0.5 * jsum * xpdf_vest[i+k*S];
            }
            else {
                sum += jsum * xpdf_vest[i+k*S];
            }
        }
        *(prmat++) = sum;
    }
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *prmat, *xpdf_vis, *xpdf_vest, *R;
#if ( ARGSCHECK==0 )
	mwSize *dims_xpdf_vis, *dims_xpdf_vest, *dims_R;
#else /* ( ARGSCHECK!=0 ) */ 
	mwSize *dims_xpdf_vis, *dims_xpdf_vest, *dims_R;
#endif /* ( ARGSCHECK!=0 ) */ 
	mwSize K, S;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<3 || nrhs>3 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_finalqtrapz:invalidNumInputs",
			"Three inputs required.");
	if ( nlhs<1 || nlhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_finalqtrapz:invalidNumOutputs",
			"One outputs required.");

	/* Get first input (XPDF_VIS, S-by-K double) */
	xpdf_vis = (double*) mxGetPr(prhs[0]);
	dims_xpdf_vis = (mwSize*) mxGetDimensions(prhs[0]);
	S = dims_xpdf_vis[0];
	K = dims_xpdf_vis[1];

	/* Get second input (XPDF_VEST, S-by-1-by-K double) */
	xpdf_vest = (double*) mxGetPr(prhs[1]);
	dims_xpdf_vest = (mwSize*) mxGetDimensions(prhs[1]);

	/* Get third input (R, 1-by-K-by-K double) */
	R = (double*) mxGetPr(prhs[2]);
	dims_R = (mwSize*) mxGetDimensions(prhs[2]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
#if ( ARGSCHECK==1 )
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:xpdf_visNotReal", "Input XPDF_VIS must be real.");
		if ( dims_xpdf_vis[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:xpdf_visWrongSize", "The first dimension of input XPDF_VIS has the wrong size (should be S).");
		if ( dims_xpdf_vis[1] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:xpdf_visWrongSize", "The second dimension of input XPDF_VIS has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:xpdf_vestNotReal", "Input XPDF_VEST must be real.");
		if ( dims_xpdf_vest[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:xpdf_vestWrongSize", "The first dimension of input XPDF_VEST has the wrong size (should be S).");
		if ( dims_xpdf_vest[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:xpdf_vestWrongSize", "The second dimension of input XPDF_VEST has the wrong size (should be 1).");
		if ( dims_xpdf_vest[2] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:xpdf_vestWrongSize", "The third dimension of input XPDF_VEST has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:RNotReal", "Input R must be real.");
		if ( dims_R[0] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:RWrongSize", "The first dimension of input R has the wrong size (should be 1).");
		if ( dims_R[1] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:RWrongSize", "The second dimension of input R has the wrong size (should be K).");
		if ( dims_R[2] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_finalqtrapz:RWrongSize", "The third dimension of input R has the wrong size (should be K).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (PRMAT, S-by-1 double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (S), (mwSize) (1), mxREAL);
	prmat = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	VestBMS_finalqtrapz(prmat, xpdf_vis, xpdf_vest, R, K, S);

}
