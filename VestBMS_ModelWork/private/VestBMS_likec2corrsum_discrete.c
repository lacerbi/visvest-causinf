#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "float.h"

/*
 * VestBMS_likec2corrsum_discrete.c
 *
 *VESTBMS_LIKEC2CORRSUM_DISCRETE Compute p(x_vis,x_vest|C=2) for DISCRETE,CORRELATED prior
 *
 * ================ INPUT VARIABLES ====================
 * PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-1] (double)
 * LIKE_VIS: visual likelihood. [S-by-K] (double)
 * LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
 *
 * ================ OUTPUT VARIABLES ==================
 * LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 26-Aug-2016 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

/* Number of discrete DISTINCT pairs of stimuli */
#define NSTIM 88

/* Summed term */
#define SUMMAND(N) (*(priorpdf2d+(N)) * *(like_vis+(N)) * *(like_vest+(N)))

void VestBMS_likec2corrsum_discrete( double *likec2, double *priorpdf2d, double *like_vis, double *like_vest, mwSize K, mwSize S )
{
    mwSize i,j,k;
    double sum,*p0,*vis0,*vest0;
    
    /* copy base pointers */    
    p0 = priorpdf2d;
    vis0 = like_vis;
    vest0 = like_vest;
    
    for (k = 0; k < K; k++) {
        for (j=0; j < K; j++) {
            priorpdf2d = p0;
            like_vis = vis0 + j*S;
            like_vest = vest0 + k*S;
            
            if ( S == NSTIM ) {
                *(likec2++) = SUMMAND(0) + SUMMAND(1) + SUMMAND(2) + SUMMAND(3) + SUMMAND(4) + SUMMAND(5) + SUMMAND(6) + SUMMAND(7) + SUMMAND(8) + SUMMAND(9) + SUMMAND(10) + SUMMAND(11) + SUMMAND(12) + SUMMAND(13) + SUMMAND(14) + SUMMAND(15) + SUMMAND(16) + SUMMAND(17) + SUMMAND(18) + SUMMAND(19) + SUMMAND(20) + SUMMAND(21) + SUMMAND(22) + SUMMAND(23) + SUMMAND(24) + SUMMAND(25) + SUMMAND(26) + SUMMAND(27) + SUMMAND(28) + SUMMAND(29) + SUMMAND(30) + SUMMAND(31) + SUMMAND(32) + SUMMAND(33) + SUMMAND(34) + SUMMAND(35) + SUMMAND(36) + SUMMAND(37) + SUMMAND(38) + SUMMAND(39) + SUMMAND(40) + SUMMAND(41) + SUMMAND(42) + SUMMAND(43) + SUMMAND(44) + SUMMAND(45) + SUMMAND(46) + SUMMAND(47) + SUMMAND(48) + SUMMAND(49) + SUMMAND(50) + SUMMAND(51) + SUMMAND(52) + SUMMAND(53) + SUMMAND(54) + SUMMAND(55) + SUMMAND(56) + SUMMAND(57) + SUMMAND(58) + SUMMAND(59) + SUMMAND(60) + SUMMAND(61) + SUMMAND(62) + SUMMAND(63) + SUMMAND(64) + SUMMAND(65) + SUMMAND(66) + SUMMAND(67) + SUMMAND(68) + SUMMAND(69) + SUMMAND(70) + SUMMAND(71) + SUMMAND(72) + SUMMAND(73) + SUMMAND(74) + SUMMAND(75) + SUMMAND(76) + SUMMAND(77) + SUMMAND(78) + SUMMAND(79) + SUMMAND(80) + SUMMAND(81) + SUMMAND(82) + SUMMAND(83) + SUMMAND(84) + SUMMAND(85) + SUMMAND(86) + SUMMAND(87);
            }
            else {
                sum = 0.;
                for (i = 0; i < S; i++)
                    sum += *(priorpdf2d++) * *(like_vis++) * *(like_vest++);
                *(likec2++) = sum;
            }
            
        }
    }
    
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *likec2, *priorpdf2d, *like_vis, *like_vest;
#if ( ARGSCHECK==0 )
	mwSize *dims_priorpdf2d, *dims_like_vis, *dims_like_vest;
#else /* ( ARGSCHECK!=0 ) */ 
	mwSize *dims_priorpdf2d, *dims_like_vis, *dims_like_vest;
#endif /* ( ARGSCHECK!=0 ) */ 
	mwSize K, S;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<3 || nrhs>3 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_likec2corrsum_discrete:invalidNumInputs",
			"Three inputs required.");
	if ( nlhs<1 || nlhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_likec2corrsum_discrete:invalidNumOutputs",
			"One outputs required.");

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

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
#if ( ARGSCHECK==1 )
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:priorpdf2dNotReal", "Input PRIORPDF2D must be real.");
		if ( dims_priorpdf2d[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:priorpdf2dWrongSize", "The first dimension of input PRIORPDF2D has the wrong size (should be S).");
		if ( dims_priorpdf2d[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:priorpdf2dWrongSize", "The second dimension of input PRIORPDF2D has the wrong size (should be 1).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:like_visNotReal", "Input LIKE_VIS must be real.");
		if ( dims_like_vis[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:like_visWrongSize", "The first dimension of input LIKE_VIS has the wrong size (should be S).");
		if ( dims_like_vis[1] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:like_visWrongSize", "The second dimension of input LIKE_VIS has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:like_vestNotReal", "Input LIKE_VEST must be real.");
		if ( dims_like_vest[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:like_vestWrongSize", "The first dimension of input LIKE_VEST has the wrong size (should be S).");
		if ( dims_like_vest[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:like_vestWrongSize", "The second dimension of input LIKE_VEST has the wrong size (should be 1).");
		if ( dims_like_vest[2] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrsum_discrete:like_vestWrongSize", "The third dimension of input LIKE_VEST has the wrong size (should be K).");

#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (LIKEC2, K-by-K double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	likec2 = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	VestBMS_likec2corrsum_discrete(likec2, priorpdf2d, like_vis, like_vest, K, S);

}
