#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "float.h"

/*
 * VestBMS_c2corrpostandlikec2adapt_mex.c
 *
 *VESTBMS_C2CORRPOSTANDLIKEC2QTRAPZ Multiple computations for C=2 (correlated) 
 *
 * ================ INPUT VARIABLES ====================
 * PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-S] (double)
 * LIKE_VIS: visual likelihood. [S-by-K] (double)
 * LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
 *
 * ================ OUTPUT VARIABLES ==================
 * POSTRIGHT_C2: p(right|x_vis,x_vest,C=2). [K-by-K] (double)
 * LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)
 * ERR: Simple estimate. [K-by-K] (double)
 * FEVALS: Function evaluations. [K-by-K] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 25-Aug-2016 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

/* Number of points per subgrid (per dimension) */
#define NGRID 5

/* Recursive adaptive integral evaluation */
double evalint(int *fevals, double *priorpdf2d, double *like_vis, double *like_vest, mwSize idx0ij, mwSize idx0i, mwSize idx0j, mwSize M, mwSize K, mwSize S, mwSize nx, mwSize ny, double tol) {
    mwSize i,j,idxi,idxj,idxij,idx0sub,Mnew;
    double tmp,y,y2,sum,sum2,maxy=-DBL_MAX,miny=DBL_MAX;
    
    /* Evaluate integrand at grid points */
    y=0.;
    idxj = idx0j + (M-1)/2;
    for (j=0; j<ny; j++) {
        idxi = idx0i + (M-1)/2;
        idxij = idx0ij + ((M-1)/2)*(S+1) + j*M*S;
        sum = 0.;
        sum2 = 0.;
        for (i=0; i<nx; i++) {
            /* Integrand */
            tmp = priorpdf2d[idxij] * like_vis[idxi] * like_vest[idxj];
            
            y += tmp;
            if ( tmp > maxy ) maxy = tmp;
            if ( tmp < miny ) miny = tmp;
            idxij += M;
            idxi += M;
        }
        /* y += sum * like_vest[idxj];
        y += sum; */
        idxj += M;
    }
    
    /* Multiply by volume element */
    y *= M*M;
    
    /* Number of function evaluations */
    /* *fevals += nx*ny; */
    
    /* Estimate if the error is larger than tolerance */
    if ( ( M > 1 ) && ((maxy - miny)*M*M > tol) ) {
        Mnew = M / (mwSize) NGRID;
        y=0.;
        for (j = 0; j<ny; j++) {
            for (i = 0; i<nx; i++) {
                /* [y(i+(j-1)*nx),fevals(i+(j-1)*nx)] = evalint(Z,idx0sub,Snew,NGRID,NGRID,col,tol); */
                y += evalint(fevals,priorpdf2d,like_vis,like_vest,idx0ij+i*M+j*S*M,idx0i+i*M,idx0j+j*M,Mnew,K,S,(mwSize) NGRID,(mwSize) NGRID,tol);
            }
        }    
    }
        
    return y;
}


/* Adaptive 2D integration */
double int2dadapt( double *err1, int *fevals1, double *priorpdf2d, double *like_vis, double *like_vest, mwSize M, mwSize K, mwSize S, mwSize nx, mwSize ny, double TolErr ) 
{
    mwSize i,j,idxi,idxj,idxij;
    double y,sum,tol;
    
    /* Compute coarse integral */
    *err1 = 0.;
    idxj = (M-1)/2;
    for (j=0; j<ny; j++) {
        idxi = (M-1)/2;
        idxij = ((M-1)/2)*(S+1) + j*M*S;
        sum = 0.;
        for (i=0; i<nx; i++) {
            sum += priorpdf2d[idxij] * like_vis[idxi];
            idxij += M;
            idxi += M;
        }        
        *err1 += sum * like_vest[idxj];
        idxj += M;
    }
    
    /* *fevals1 += (int) (nx*ny); */    
    *err1 *= (double) (M*M);
    tol = (*err1) * TolErr;    /* relative error tolerance */
    
    y = evalint(fevals1,priorpdf2d,like_vis,like_vest,0,0,0,M,K,S,nx,ny,tol);

    return y;
}




void VestBMS_c2corrpostandlikec2adapt( double *postright_c2, double *likec2, double *err, double *fevals, double *priorpdf2d, double *like_vis, double *like_vest, mwSize K, mwSize S )
{
    mwSize M,nx,ny,k,m;
    double pmin,TolErr,y1,y2,err1,err2;
    int fevals1,fevals2;

    /* Accepted relative error (approximate) */
    TolErr = 0.001;
    
    M = (mwSize) (NGRID*NGRID*NGRID);
    nx = S / M;
    ny = S / M;
    
    /* Add some small probability to prevent 0 / 0 */
    pmin = 0.5 * DBL_MIN * (double) S;    

    /* Outer loop */
    for (m = 0; m < K; m++) {
        for (k = 0; k < K; k++) {

            fevals1=0;
            fevals2=0;
            
            y1 = int2dadapt(&err1,&fevals1,priorpdf2d,&like_vis[k*S],&like_vest[m*S],M,K,S,nx,ny/2,TolErr);
            y2 = int2dadapt(&err2,&fevals2,&priorpdf2d[S*S/2],&like_vis[k*S],&like_vest[m*S+(S/2)],M,K,S,nx,ny/2,TolErr);
            
            err[k + m*K] = err1 + err2;
            fevals[k + m*K] = (double) (fevals1 + fevals2);

            likec2[k + m*K] = y1 + y2 + pmin;
            postright_c2[k + m*K] = (y2 + pmin) / (y1 + y2 + pmin);

        }
    }
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *postright_c2, *likec2, *err, *priorpdf2d, *like_vis, *like_vest;
	int *fevals;
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
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_c2corrpostandlikec2adapt:invalidNumInputs",
			"Three inputs required.");
	if ( nlhs<4 || nlhs>4 )
		mexErrMsgIdAndTxt( "MATLAB:VestBMS_c2corrpostandlikec2adapt:invalidNumOutputs",
			"Four outputs required.");

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
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:priorpdf2dNotReal", "Input PRIORPDF2D must be real.");
		if ( dims_priorpdf2d[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:priorpdf2dWrongSize", "The first dimension of input PRIORPDF2D has the wrong size (should be S).");
		if ( dims_priorpdf2d[1] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:priorpdf2dWrongSize", "The second dimension of input PRIORPDF2D has the wrong size (should be S).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:like_visNotReal", "Input LIKE_VIS must be real.");
		if ( dims_like_vis[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:like_visWrongSize", "The first dimension of input LIKE_VIS has the wrong size (should be S).");
		if ( dims_like_vis[1] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:like_visWrongSize", "The second dimension of input LIKE_VIS has the wrong size (should be K).");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:like_vestNotReal", "Input LIKE_VEST must be real.");
		if ( dims_like_vest[0] != ((mwSize) (S)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:like_vestWrongSize", "The first dimension of input LIKE_VEST has the wrong size (should be S).");
		if ( dims_like_vest[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:like_vestWrongSize", "The second dimension of input LIKE_VEST has the wrong size (should be 1).");
		if ( dims_like_vest[2] != ((mwSize) (K)) )
			mexErrMsgIdAndTxt("MATLAB:VestBMS_c2corrpostandlikec2adapt:like_vestWrongSize", "The third dimension of input LIKE_VEST has the wrong size (should be K).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (POSTRIGHT_C2, K-by-K double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	postright_c2 = mxGetPr(plhs[0]);

	/* Pointer to second output (LIKEC2, K-by-K double) */
	plhs[1] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	likec2 = mxGetPr(plhs[1]);

	/* Pointer to third output (ERR, K-by-K double) */
	plhs[2] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	err = mxGetPr(plhs[2]);

	/* Pointer to fourth output (FEVALS, K-by-K int) */
	plhs[3] = mxCreateDoubleMatrix((mwSize) (K), (mwSize) (K), mxREAL);
	fevals = mxGetPr(plhs[3]);

	/* Call the C subroutine */
	VestBMS_c2corrpostandlikec2adapt(postright_c2, likec2, err, fevals, priorpdf2d, like_vis, like_vest, K, S);

}
