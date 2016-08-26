#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "float.h"

/*
 * VestBMS_c1postandlikec1sum_discrete.c
 *
 *VESTBMS_C1POSTANDLIKEC1QTRAPZ Multiple computations for C=1 (correlated,discrete)
 *
 * ================ INPUT VARIABLES ====================
 * POSTPDF_C2_UNI: p(s) * p(x_vest|s). [S-by-1-by-K] (double)
 * LIKE_VIS_UNI: p(x_vis|s). [S-by-K] (double)
 * SRANGE_VEST_UNI: s range. [S-by-1] (double)
 *
 * ================ OUTPUT VARIABLES ==================
 * POSTRIGHT_C1: p(right|x_vis,x_vest,C=1). [K-by-K] (double)
 * LIKEC1: p(x_vis,x_vest|C=1). [K-by-K] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 26-Aug-2016 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void VestBMS_c1postandlikec1sum_discrete( double *postright_c1, double *likec1, double *postpdf_c2_uni, double *like_vis_uni, double *srange_uni, mwSize K, mwSize S )
{
    mwSize i,j,k;
    double pmin,sum0,suml,sumr,*p0,*vis0,*vest0;

    /* Add some small probability to prevent 0 / 0 */
    pmin = 0.5 * DBL_MIN * (double) S;    
    
    /* copy base pointers */    
    p0 = postpdf_c2_uni;
    vis0 = like_vis_uni;
    
    for (j=0; j < K; j++) {
        for (k = 0; k < K; k++) {
            postpdf_c2_uni = p0 + k*S;
            like_vis_uni = vis0 + j*S;
            sum0 = 0.;
            suml = 0.;
            sumr = 0.;
            
            for (i = 0; i < S; i++) {
                if ( srange_uni[i] < 0. ) {
                    suml += *(postpdf_c2_uni++) * *(like_vis_uni++);
                }
                else if ( srange_uni[i] > 0. ) {
                    sumr += *(postpdf_c2_uni++) * *(like_vis_uni++);
                }
                else {
                    sum0 += *(postpdf_c2_uni++) * *(like_vis_uni++);
                }
            }
            
            likec1[k*K+j] = suml + sumr + sum0 + pmin;
            postright_c1[k*K+j] = (sumr + 0.5*sum0 + pmin) / likec1[k*K+j];
            
        }
    }
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *postright_c1, *likec1, *postpdf_c2_uni, *like_vis_uni, *srange_uni;
#if ( ARGSCHECK==0 )
	mwSize *dims_postpdf_c2_uni, *dims_like_vis_uni, *dims_srange_uni;
#else /* ( ARGSCHECK!=0 ) */ 
	mwSize *dims_postpdf_c2_uni, *dims_like_vis_uni, *dims_srange_uni;
#endif /* ( ARGSCHECK!=0 ) */ 
	mwSize K, S;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<3 || nrhs>3 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_c1postandlikec1sum_discrete:invalidNumInputs",
			"Three inputs required.");
	if ( nlhs<2 || nlhs>2 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_c1postandlikec1sum_discrete:invalidNumOutputs",
			"Two outputs required.");

	/* Get first input (POSTPDF_C2_UNI, S-by-1-by-K double) */
	postpdf_c2_uni = (double*) mxGetPr(prhs[0]);
	dims_postpdf_c2_uni = (mwSize*) mxGetDimensions(prhs[0]);
	S = dims_postpdf_c2_uni[0];
	K = dims_postpdf_c2_uni[2];

	/* Get second input (LIKE_VIS_UNI, S-by-K double) */
	like_vis_uni = (double*) mxGetPr(prhs[1]);
	dims_like_vis_uni = (mwSize*) mxGetDimensions(prhs[1]);

	/* Get third input (SRANGE_UNI, S-by-1 double) */
	srange_uni = (double*) mxGetPr(prhs[2]);
	dims_srange_uni = (mwSize*) mxGetDimensions(prhs[2]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
#if ( ARGSCHECK==1 )
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:postpdf_c2_uniNotReal", "Input POSTPDF_C2_UNI must be real.");
		if ( dims_postpdf_c2_uni[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:postpdf_c2_uniWrongSize", "The first dimension of input POSTPDF_C2_UNI has the wrong size (should be S).");
		if ( dims_postpdf_c2_uni[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:postpdf_c2_uniWrongSize", "The second dimension of input POSTPDF_C2_UNI has the wrong size (should be 1).");
		if ( dims_postpdf_c2_uni[2] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:postpdf_c2_uniWrongSize", "The third dimension of input POSTPDF_C2_UNI has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:like_vis_uniNotReal", "Input LIKE_VIS_UNI must be real.");
		if ( dims_like_vis_uni[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:like_vis_uniWrongSize", "The first dimension of input LIKE_VIS_UNI has the wrong size (should be S).");
		if ( dims_like_vis_uni[1] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:like_vis_uniWrongSize", "The second dimension of input LIKE_VIS_UNI has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:srange_uniNotReal", "Input SRANGE_UNI must be real.");
		if ( dims_srange_uni[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:srange_uniWrongSize", "The first dimension of input SRANGE_UNI has the wrong size (should be S).");
		if ( dims_srange_uni[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c1postandlikec1sum_discrete:srange_uniWrongSize", "The second dimension of input SRANGE_UNI has the wrong size (should be 1).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (POSTRIGHT_C1, K-by-K double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	postright_c1 = mxGetPr(plhs[0]);

	/* Pointer to second output (LIKEC1, K-by-K double) */
	plhs[1] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	likec1 = mxGetPr(plhs[1]);

	/* Call the C subroutine */
	VestBMS_c1postandlikec1sum_discrete(postright_c1, likec1, postpdf_c2_uni, like_vis_uni, srange_uni, K, S);

}
