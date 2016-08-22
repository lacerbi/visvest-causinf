#include "mex.h"
#include "math.h"
#include "matrix.h"

/*
 * VestBMS_likec2corrqtrapz_mat_mex.c
 *
 *VESTBMS_LIKEC2CORRQTRAPZ Compute p(x_vis,x_vest|C=2) for CORRELATED prior
 *
 * ================ INPUT VARIABLES ====================
 * PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-S] (double)
 * LIKE_VIS: visual likelihood. [S-by-K] (double)
 * LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
 *
 * ================ OUTPUT VARIABLES ==================
 * LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 21-Aug-2016 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void VestBMS_likec2corrqtrapz( double *likec2, double *priorpdf2d, double *like_vis, double *like_vest, mwSize K, mwSize S )
{
	
    mwSize i,j,k,l;
    double SUM,sum,*p0,*vis0,*vest0;
    
    /* copy base pointers */    
    p0 = priorpdf2d;
    vis0 = like_vis;
    vest0 = like_vest;
    
    for (k=0; k<K; k++) {
        for (l=0; l<K; l++) {
            like_vest = vest0 + l*S;            
            priorpdf2d = p0;

            /* case j = 0 */
            like_vis = vis0 + k*S;
            sum = 0.25 * *(priorpdf2d++) * *(like_vis++);
            for (i=1; i<S-1; i++) {
                sum += 0.5 * *(priorpdf2d++) * *(like_vis++);
            }
            sum += 0.25 * *(priorpdf2d++) * *(like_vis++);
            SUM = sum * *(like_vest++);     /* Initialize SUM */
            
            /* case j = 1 to S-2 */            
            for (j=1; j<S-1; j++) {
                like_vis = vis0 + k*S;
                
                sum = 0.5 * *(priorpdf2d++) * *(like_vis++);
                for (i=1; i<S-1; i++) {
                    sum += *(priorpdf2d++) * *(like_vis++);
                }
                sum += 0.5 * *(priorpdf2d++) * *(like_vis++);
                SUM += sum * *(like_vest++);
            }

            /* case j = S-1 */
            like_vis = vis0 + k*S;
            sum = 0.25 * *(priorpdf2d++) * *(like_vis++);
            for (i=1; i<S-1; i++) {
                sum += 0.5 * *(priorpdf2d++) * *(like_vis++);
            }
            sum += 0.25 * *(priorpdf2d++) * *(like_vis++);
            SUM += sum * *(like_vest++);
            
            likec2[l*K+k] = SUM;

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
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_likec2corrqtrapz:invalidNumInputs",
			"Three inputs required.");
	if ( nlhs<1 || nlhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_likec2corrqtrapz:invalidNumOutputs",
			"One outputs required.");

	/* Get first input (PRIORPDF2D, S-by-S double) */
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
				mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:priorpdf2dNotReal", "Input PRIORPDF2D must be real.");
		if ( dims_priorpdf2d[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:priorpdf2dWrongSize", "The first dimension of input PRIORPDF2D has the wrong size (should be S).");
		if ( dims_priorpdf2d[1] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:priorpdf2dWrongSize", "The second dimension of input PRIORPDF2D has the wrong size (should be S).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:like_visNotReal", "Input LIKE_VIS must be real.");
		if ( dims_like_vis[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:like_visWrongSize", "The first dimension of input LIKE_VIS has the wrong size (should be S).");
		if ( dims_like_vis[1] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:like_visWrongSize", "The second dimension of input LIKE_VIS has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:like_vestNotReal", "Input LIKE_VEST must be real.");
		if ( dims_like_vest[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:like_vestWrongSize", "The first dimension of input LIKE_VEST has the wrong size (should be S).");
		if ( dims_like_vest[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:like_vestWrongSize", "The second dimension of input LIKE_VEST has the wrong size (should be 1).");
		if ( dims_like_vest[2] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_likec2corrqtrapz:like_vestWrongSize", "The third dimension of input LIKE_VEST has the wrong size (should be K).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (LIKEC2, K-by-K double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	likec2 = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	VestBMS_likec2corrqtrapz(likec2, priorpdf2d, like_vis, like_vest, K, S);

}
