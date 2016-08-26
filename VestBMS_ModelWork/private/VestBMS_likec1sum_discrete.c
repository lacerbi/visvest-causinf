#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "float.h"

/*
 * VestBMS_likec1sum_discrete.c
 *
 *VESTBMS_LIKEC1SUM_DISCRETE Multiple computations for C=1 (correlated,discrete)
 *
 * ================ INPUT VARIABLES ====================
 * POSTPDF_C2_UNI: p(s) * p(x_vest|s). [S-by-1-by-K] (double)
 * LIKE_VIS_UNI: p(x_vis|s). [S-by-K] (double)
 *
 * ================ OUTPUT VARIABLES ==================
 * LIKEC1: p(x_vis,x_vest|C=1). [K-by-K] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 26-Aug-2016 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void VestBMS_likec1sum_discrete( double *likec1, double *postpdf_c2_uni, double *like_vis_uni, mwSize K, mwSize S )
{
    mwSize i,j,k;
    double pmin,sum,*p0,*vis0;
    
    /* copy base pointers */    
    p0 = postpdf_c2_uni;
    vis0 = like_vis_uni;
    
    for (j=0; j < K; j++) {
        for (k = 0; k < K; k++) {
            postpdf_c2_uni = p0 + k*S;
            like_vis_uni = vis0 + j*S;
            sum = 0.;            
            for (i = 0; i < S; i++)
                sum += *(postpdf_c2_uni++) * *(like_vis_uni++);            
            likec1[k*K+j] = sum;            
        }
    }
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *likec1, *postpdf_c2_uni, *like_vis_uni;
#if ( ARGSCHECK==0 )
	mwSize *dims_postpdf_c2_uni, *dims_like_vis_uni;
#else /* ( ARGSCHECK!=0 ) */ 
	mwSize *dims_postpdf_c2_uni, *dims_like_vis_uni;
#endif /* ( ARGSCHECK!=0 ) */ 
	mwSize K, S;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<2 || nrhs>2 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_likec1sum_discrete:invalidNumInputs",
			"Two inputs required.");
	if ( nlhs<1 || nlhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_likec1sum_discrete:invalidNumOutputs",
			"One outputs required.");

	/* Get first input (POSTPDF_C2_UNI, S-by-1-by-K double) */
	postpdf_c2_uni = (double*) mxGetPr(prhs[0]);
	dims_postpdf_c2_uni = (mwSize*) mxGetDimensions(prhs[0]);
	S = dims_postpdf_c2_uni[0];
	K = dims_postpdf_c2_uni[2];

	/* Get second input (LIKE_VIS_UNI, S-by-K double) */
	like_vis_uni = (double*) mxGetPr(prhs[1]);
	dims_like_vis_uni = (mwSize*) mxGetDimensions(prhs[1]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
#if ( ARGSCHECK==1 )
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_likec1sum_discrete:postpdf_c2_uniNotReal", "Input POSTPDF_C2_UNI must be real.");
		if ( dims_postpdf_c2_uni[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec1sum_discrete:postpdf_c2_uniWrongSize", "The first dimension of input POSTPDF_C2_UNI has the wrong size (should be S).");
		if ( dims_postpdf_c2_uni[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec1sum_discrete:postpdf_c2_uniWrongSize", "The second dimension of input POSTPDF_C2_UNI has the wrong size (should be 1).");
		if ( dims_postpdf_c2_uni[2] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec1sum_discrete:postpdf_c2_uniWrongSize", "The third dimension of input POSTPDF_C2_UNI has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_likec1sum_discrete:like_vis_uniNotReal", "Input LIKE_VIS_UNI must be real.");
		if ( dims_like_vis_uni[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec1sum_discrete:like_vis_uniWrongSize", "The first dimension of input LIKE_VIS_UNI has the wrong size (should be S).");
		if ( dims_like_vis_uni[1] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec1sum_discrete:like_vis_uniWrongSize", "The second dimension of input LIKE_VIS_UNI has the wrong size (should be K).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (LIKEC1, K-by-K double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	likec1 = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	VestBMS_likec1sum_discrete(likec1, postpdf_c2_uni, like_vis_uni, K, S);

}
