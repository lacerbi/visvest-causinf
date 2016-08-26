#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "float.h"

/*
 * VestBMS_c2corrpostandlikec2sum_discrete_mat_mex.c
 *
 *VESTBMS_C2CORRPOSTANDLIKEC2SUM_DISCRETE Multiple computations for C=2 (correlated,discrete)
 *
 * ================ INPUT VARIABLES ====================
 * PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-1] (double)
 * LIKE_VIS: visual likelihood. [S-by-K] (double)
 * LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
 * SRANGE_VEST: vestibular stimuli. [S-by-1] (double)
 *
 * ================ OUTPUT VARIABLES ==================
 * POSTRIGHT_C2: p(right|x_vis,x_vest,C=2). [K-by-K] (double)
 * LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 26-Aug-2016 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void VestBMS_c2corrpostandlikec2sum_discrete( double *postright_c2, double *likec2, double *priorpdf2d, double *like_vis, double *like_vest, double *srange_vest, mwSize K, mwSize S )
{
    mwSize i,j,k;
    double pmin,sum0,suml,sumr,*p0,*vis0,*vest0;

    /* Add some small probability to prevent 0 / 0 */
    pmin = 0.5 * DBL_MIN * (double) S;    
    
    /* copy base pointers */    
    p0 = priorpdf2d;
    vis0 = like_vis;
    vest0 = like_vest;
    
    for (j=0; j < K; j++) {
        for (k = 0; k < K; k++) {
            priorpdf2d = p0;
            like_vis = vis0 + j*S;
            like_vest = vest0 + k*S;
            sum0 = 0.;
            suml = 0.;
            sumr = 0.;
            
            for (i = 0; i < S; i++) {
                if ( srange_vest[i] < 0. ) {
                    suml += *(priorpdf2d++) * *(like_vis++) * *(like_vest++);
                }
                else if ( srange_vest[i] > 0. ) {
                    sumr += *(priorpdf2d++) * *(like_vis++) * *(like_vest++);
                }
                else {
                    sum0 += *(priorpdf2d++) * *(like_vis++) * *(like_vest++);
                }
            }
            
            likec2[k*K+j] = suml + sumr + sum0 + pmin;
            postright_c2[k*K+j] = (sumr + 0.5*sum0 + pmin) / likec2[k*K+j];
            
        }
    }
    
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *postright_c2, *likec2, *priorpdf2d, *like_vis, *like_vest, *srange_vest;
#if ( ARGSCHECK==0 )
	mwSize *dims_priorpdf2d, *dims_like_vis, *dims_like_vest, *dims_srange_vest;
#else /* ( ARGSCHECK!=0 ) */ 
	mwSize *dims_priorpdf2d, *dims_like_vis, *dims_like_vest, *dims_srange_vest;
#endif /* ( ARGSCHECK!=0 ) */ 
	mwSize K, S;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<4 || nrhs>4 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:invalidNumInputs",
			"Four inputs required.");
	if ( nlhs<2 || nlhs>2 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:invalidNumOutputs",
			"Two outputs required.");

	/* Get first input (PRIORPDF2D, S-by-1 double) */
	priorpdf2d = (double*) mxGetPr(prhs[0]);
	dims_priorpdf2d = (mwSize*) mxGetDimensions(prhs[0]);
	S = dims_priorpdf2d[0];

	/* Get second input (LIKE_VIS, S-by-K double) */
	like_vis = (double*) mxGetPr(prhs[1]);
	dims_like_vis = (mwSize*) mxGetDimensions(prhs[1]);
	K = dims_like_vis[1];

	/* Get third input (LIKE_VEST, S-by-1-by-K double) */
	like_vest = (double*) mxGetPr(prhs[2]);
	dims_like_vest = (mwSize*) mxGetDimensions(prhs[2]);

	/* Get fourth input (SRANGE_VEST, S-by-1 double) */
	srange_vest = (double*) mxGetPr(prhs[3]);
	dims_srange_vest = (mwSize*) mxGetDimensions(prhs[3]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
#if ( ARGSCHECK==1 )
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:priorpdf2dNotReal", "Input PRIORPDF2D must be real.");
		if ( dims_priorpdf2d[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:priorpdf2dWrongSize", "The first dimension of input PRIORPDF2D has the wrong size (should be S).");
		if ( dims_priorpdf2d[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:priorpdf2dWrongSize", "The second dimension of input PRIORPDF2D has the wrong size (should be 1).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:like_visNotReal", "Input LIKE_VIS must be real.");
		if ( dims_like_vis[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:like_visWrongSize", "The first dimension of input LIKE_VIS has the wrong size (should be S).");
		if ( dims_like_vis[1] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:like_visWrongSize", "The second dimension of input LIKE_VIS has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:like_vestNotReal", "Input LIKE_VEST must be real.");
		if ( dims_like_vest[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:like_vestWrongSize", "The first dimension of input LIKE_VEST has the wrong size (should be S).");
		if ( dims_like_vest[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:like_vestWrongSize", "The second dimension of input LIKE_VEST has the wrong size (should be 1).");
		if ( dims_like_vest[2] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:like_vestWrongSize", "The third dimension of input LIKE_VEST has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:srange_vestNotReal", "Input SRANGE_VEST must be real.");
		if ( dims_srange_vest[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:srange_vestWrongSize", "The first dimension of input SRANGE_VEST has the wrong size (should be S).");
		if ( dims_srange_vest[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2sum_discrete:srange_vestWrongSize", "The second dimension of input SRANGE_VEST has the wrong size (should be 1).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (POSTRIGHT_C2, K-by-K double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	postright_c2 = mxGetPr(plhs[0]);

	/* Pointer to second output (LIKEC2, K-by-K double) */
	plhs[1] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	likec2 = mxGetPr(plhs[1]);

	/* Call the C subroutine */
	VestBMS_c2corrpostandlikec2sum_discrete(postright_c2, likec2, priorpdf2d, like_vis, like_vest, srange_vest, K, S);

}
