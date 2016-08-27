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

/* Number of discrete pairs of DISTINCT stimuli */
#define NSTIM 88

/* Summed term */
#define SUMMAND(N) (*(priorpdf2d+(N)) * *(like_vis+(N)) * *(like_vest+(N)))

/* Hardcoded version much faster */
void VestBMS_c2corrpostandlikec2sum_discrete_hardcoded( double *postright_c2, double *likec2, double *priorpdf2d, double *like_vis, double *like_vest, mwSize K, mwSize S )
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
            
            suml = SUMMAND(0) + SUMMAND(11) + SUMMAND(12) + SUMMAND(13) + SUMMAND(22) + SUMMAND(23) + SUMMAND(24) + SUMMAND(25) + SUMMAND(33) + SUMMAND(34) + SUMMAND(35) + SUMMAND(36) + SUMMAND(37) + SUMMAND(44) + SUMMAND(45) + SUMMAND(46) + SUMMAND(47) + SUMMAND(48) + SUMMAND(49) + SUMMAND(55) + SUMMAND(56) + SUMMAND(57) + SUMMAND(58) + SUMMAND(59) + SUMMAND(60) + SUMMAND(66) + SUMMAND(67) + SUMMAND(68) + SUMMAND(69) + SUMMAND(70) + SUMMAND(71) + SUMMAND(72) + SUMMAND(77) + SUMMAND(78) + SUMMAND(79) + SUMMAND(80) + SUMMAND(81) + SUMMAND(82) + SUMMAND(83) + SUMMAND(84) + SUMMAND(85);
            sumr = SUMMAND(2) + SUMMAND(3) + SUMMAND(4) + SUMMAND(5) + SUMMAND(6) + SUMMAND(7) + SUMMAND(8) + SUMMAND(9) + SUMMAND(10) + SUMMAND(15) + SUMMAND(16) + SUMMAND(17) + SUMMAND(18) + SUMMAND(19) + SUMMAND(20) + SUMMAND(21) + SUMMAND(27) + SUMMAND(28) + SUMMAND(29) + SUMMAND(30) + SUMMAND(31) + SUMMAND(32) + SUMMAND(38) + SUMMAND(39) + SUMMAND(40) + SUMMAND(41) + SUMMAND(42) + SUMMAND(43) + SUMMAND(50) + SUMMAND(51) + SUMMAND(52) + SUMMAND(53) + SUMMAND(54) + SUMMAND(62) + SUMMAND(63) + SUMMAND(64) + SUMMAND(65) + SUMMAND(74) + SUMMAND(75) + SUMMAND(76) + SUMMAND(87);
            sum0 = SUMMAND(1) + SUMMAND(14) + SUMMAND(26) + SUMMAND(61) + SUMMAND(73) + SUMMAND(86);
            
            /* This sum incorrectly contains terms with C=1 
             suml = SUMMAND(0) + SUMMAND(11) + SUMMAND(12) + SUMMAND(13) + SUMMAND(22) + SUMMAND(23) + SUMMAND(24) + SUMMAND(25) + SUMMAND(33) + SUMMAND(34) + SUMMAND(35) + SUMMAND(36) + SUMMAND(37) + SUMMAND(44) + SUMMAND(45) + SUMMAND(46) + SUMMAND(47) + SUMMAND(48) + SUMMAND(55) + SUMMAND(56) + SUMMAND(57) + SUMMAND(58) + SUMMAND(59) + SUMMAND(60) + SUMMAND(66) + SUMMAND(67) + SUMMAND(68) + SUMMAND(69) + SUMMAND(70) + SUMMAND(71) + SUMMAND(77) + SUMMAND(78) + SUMMAND(79) + SUMMAND(80) + SUMMAND(81) + SUMMAND(82) + SUMMAND(83) + SUMMAND(88) + SUMMAND(89) + SUMMAND(90) + SUMMAND(91) + SUMMAND(92) + SUMMAND(93) + SUMMAND(94) + SUMMAND(95) + SUMMAND(96);
            sumr = SUMMAND(2) + SUMMAND(3) + SUMMAND(4) + SUMMAND(5) + SUMMAND(6) + SUMMAND(7) + SUMMAND(8) + SUMMAND(9) + SUMMAND(10) + SUMMAND(15) + SUMMAND(16) + SUMMAND(17) + SUMMAND(18) + SUMMAND(19) + SUMMAND(20) + SUMMAND(21) + SUMMAND(27) + SUMMAND(28) + SUMMAND(29) + SUMMAND(30) + SUMMAND(31) + SUMMAND(32) + SUMMAND(38) + SUMMAND(39) + SUMMAND(40) + SUMMAND(41) + SUMMAND(42) + SUMMAND(43) + SUMMAND(50) + SUMMAND(51) + SUMMAND(52) + SUMMAND(53) + SUMMAND(54) + SUMMAND(61) + SUMMAND(62) + SUMMAND(63) + SUMMAND(64) + SUMMAND(65) + SUMMAND(73) + SUMMAND(74) + SUMMAND(75) + SUMMAND(76) + SUMMAND(85) + SUMMAND(86) + SUMMAND(87) + SUMMAND(98);
            sum0 = SUMMAND(1) + SUMMAND(14) + SUMMAND(26) + SUMMAND(49) + SUMMAND(72) + SUMMAND(84) + SUMMAND(97); */
            
            likec2[k*K+j] = suml + sumr + sum0 + pmin;
            postright_c2[k*K+j] = (sumr + 0.5*sum0 + pmin) / likec2[k*K+j];            
        }
    }
    
	
}

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
	mwSize i, K, S;
    /* Hardcoded stimulus vector */
    double srange[NSTIM] = {-5.0,0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,-15.0,-10.0,-5.0,0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,-20.0,-15.0,-10.0,-5.0,0.0,5.0,10.0,15.0,20.0,25.0,30.0,-22.5,-17.5,-12.5,-7.5,-2.5,2.5,7.5,12.5,17.5,22.5,27.5,-27.5,-22.5,-17.5,-12.5,-7.5,-2.5,2.5,7.5,12.5,17.5,22.5,-30.0,-25.0,-20.0,-15.0,-10.0,-5.0,0.0,5.0,10.0,15.0,20.0,-35.0,-30.0,-25.0,-20.0,-15.0,-10.0,-5.0,0.0,5.0,10.0,15.0,-45.0,-40.0,-35.0,-30.0,-25.0,-20.0,-15.0,-10.0,-5.0,0.0,5.0};

    /* This vector erroneously contained also the C=1 cases */
    /* double srange[NSTIM] = {-5.,0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.,35.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.,-22.5,-17.5,-12.5,-7.5,-2.5,2.5,7.5,12.5,17.5,22.5,27.5,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,-27.5,-22.5,-17.5,-12.5,-7.5,-2.5,2.5,7.5,12.5,17.5,22.5,-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,-35.,-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,-45.,-40.,-35.,-30.,-25.,-20.,-15.,-10.,-5.,0.,5.}; */

    int hard;

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

    /* Check whether to call normal function or hard-coded */
    hard = 1;
    if ( S != NSTIM ) hard = 0;
    else {    
        for (i=0; i<99; i++) {
            if (srange[i] != srange_vest[i]) {
                hard = 0;
                break;
            }
        }
    }
        
	/* Call the C subroutine */
    if ( hard )
        VestBMS_c2corrpostandlikec2sum_discrete_hardcoded(postright_c2, likec2, priorpdf2d, like_vis, like_vest, K, S);
    else
        VestBMS_c2corrpostandlikec2sum_discrete(postright_c2, likec2, priorpdf2d, like_vis, like_vest, srange_vest, K, S);
}
