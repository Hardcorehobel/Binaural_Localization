//  xcorrNormMEX   Normalized time-domain cross-correlation function.
//
//  Author  :  Tobias May, © 2007-2009 
//             TUe Eindhoven and Philips Research  
//             t.may@tue.nl      tobias.may@philips.com
//
//  History :  
//  v.0.1   2007/11/08
//  v.0.2   2009/10/11 cleaned up
//  ***********************************************************************

#include <stdio.h>
#include "math.h"
#include "mex.h"

/* Input Arguments */
#define   INPUT1      			 prhs[0]  // left input signal 
#define   INPUT2      			 prhs[1]  // right input signal 
#define   MAXLAG    			 prhs[2]  // maximum lag
#define   DETREND                prhs[3]  // detrend flag
#define   NORM                   prhs[4]  // normalization flag

/* Output Arguments */
#define   OUTPUT				 plhs[0]  // output signal
#define   LAGS  				 plhs[1]  // lags

/* Helper functions */
#define max(x, y)   ((x) > (y) ? (x) : (y))
#define	min(A, B)	((A) < (B) ? (A) : (B))
#define swap(A,B)   temp = (A); (A)=(B); (B) = temp;
#define getRound(x) ((x) >= 0?(long)((x)+0.5):(long)((x)-0.5))

/*  Set resolution of floating point numbers (equivalent to MATLABs eps) */
const double EPS = pow(2.00,-52.00);

// Help 
void usage()
{
	mexPrintf(" xcorrNorm   Normalized time-domain cross-correlation function.\n"); 
	mexPrintf("\n"); 
	mexPrintf(" USAGE\n"); 
	mexPrintf("\t [XCORR,LAGS] = xcorrNorm(INL,INR)\n"); 
	mexPrintf("\t [XCORR,LAGS] = xcorrNorm(INL,INR,MAXLAG,bDETREND,bNORM)\n"); 
	mexPrintf("\n"); 
	mexPrintf(" INPUT ARGUMENTS\n");
	mexPrintf("\t      INL : left input arranged as  [nSamples x nChannels]\n");
	mexPrintf("\t      INR : right input arranged as [nSamples x nChannels]\n");
	mexPrintf("\t   MAXLAG : computation is performned over the lag range -MAXLAG:MAXLAG\n");
	mexPrintf("\t            (default, MAXLAG = nSamples-1) \n");
	mexPrintf("\t bDETREND : substract mean     (default, bDETREND = true) \n");
	mexPrintf("\t    bNORM : normalization flag (default, bNORM    = true) \n");
	mexPrintf("\n");
	mexPrintf(" OUTPUT ARGUMENTS\n");
	mexPrintf("\t   XCORR : cross-correlation function [nSamples x nChannels]  \n");
	mexPrintf("\t    LAGS : time lags of cross-correlation function [2*MAXLAG+1 x 1] \n");
	mexPrintf("\n");
	mexPrintf(" REFERENCES\n");
	mexPrintf("\t [1]  Roman, N., Wang, D. L. and Brown, G. J., \"Speech segregation based\n");
	mexPrintf("\t      on sound localization\", J. Acoust. Soc. Amer., vol. 114, no. 4,\n");
	mexPrintf("\t      pp. 2236-2252, 2003.\n");
	mexPrintf("\n");
	mexPrintf("\t Author  :  Tobias May, © 2007-2009 \n");
	mexPrintf("\t            TUe Eindhoven and Philips Research   \n");
	mexPrintf("\t            t.may@tue.nl      tobias.may@philips.com \n");
	mexPrintf("\n"); 
}

// MEX wrapper
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *input1, *input2, *output, *lags;

	int    nSamples, nChannels;
	int    chan, bNorm, bDetrend, maxLag, nLags;

	double d1, d2, avg1, avg2 , cc;
	double *delayedLeft, *delayedRight;
	int    delay, win;

	// Check for proper number of arguments
	if (nrhs < 2){
		usage(); mexErrMsgTxt("Not enough input arguments.");
	}

	if (nrhs > 5){
		usage(); mexErrMsgTxt("Too many input arguments.");
	}

	if (nlhs > 2){
		usage(); mexErrMsgTxt("Too many output arguments.");
	}

	// Check if dimensions of both input signals match
	if(((int) mxGetM(INPUT1) != (int) mxGetM(INPUT2)) || ((int) mxGetN(INPUT1) != (int) mxGetN(INPUT2)) ){
		usage(); mexErrMsgTxt("Dimension mismatch between \"INL\" and \"INR\".");
	}

	// Get dimension of input data
	nSamples  = (int)mxGetM(INPUT1);
	nChannels = (int)mxGetN(INPUT1);

	// Set defaults values
	if (nrhs < 5)
		bNorm = true;
	else
		bNorm = (int)mxGetScalar(NORM);

	if (nrhs < 4)
		bDetrend = true;
	else
		bDetrend = (int)mxGetScalar(DETREND);

	if (nrhs < 3)
		maxLag = nSamples - 1;
	else
		// Limit the maximum lag to be at least 1
		maxLag = max(1, (int)mxGetScalar(MAXLAG));

	nLags = 2 * maxLag + 1;

	// Asign pointers
	input1 = mxGetPr(INPUT1);
	input2 = mxGetPr(INPUT2);

	// Create a matrix for the return argument
	OUTPUT = mxCreateDoubleMatrix(nLags, nChannels, mxREAL);

	// Asign pointers
	output = mxGetPr(OUTPUT);

	// Generate vector of lags
	if (nlhs > 1){
		LAGS = mxCreateDoubleMatrix(nLags, 1, mxREAL);
		lags = mxGetPr(LAGS);
		for(int ii = 0; ii < nLags; ii++){
			lags[ii] = (double) -maxLag + ii;
		}
	}

	delayedLeft  = (double*) mxCalloc(nSamples, sizeof(double));
	delayedRight = (double*) mxCalloc(nSamples, sizeof(double));

	// Loop over number of frames
	for (chan = 0; chan < nChannels; chan++){

		// Loop over number of delays
		for (delay = -maxLag; delay <=maxLag; delay++){

			avg1 = 0.0; // Mean left
			avg2 = 0.0; // Mean right
			d1   = 0.0; // Auto-correlation left
			d2   = 0.0; // Auto-correlation right
			cc   = 0.0; // Cross-product


			// ===========================================================
			// Built delayed signal vectors
			// ===========================================================
			if (delay < 0){
				for (win=0;win<nSamples;win++){
					// Shift to negative time lags
					if(win>=abs(delay))
						delayedLeft[win] = input1[chan * nSamples + win-abs(delay)];
					else
						delayedLeft[win] = 0.0;

					delayedRight[win] = input2[chan * nSamples + win];
				}
			}

			else if(delay == 0){
				for (win=0;win<nSamples;win++){
					// No time delay
					delayedLeft[win]  = input1[chan * nSamples + win];
					delayedRight[win] = input2[chan * nSamples + win];
				}
			}

			else{
				for (win=0;win<nSamples;win++){
					// Shift to positibe time lags
					if(win>=abs(delay))
						delayedRight[win] = input2[chan * nSamples + win-abs(delay)];
					else
						delayedRight[win] = 0.0;

					delayedLeft[win] = input1[chan * nSamples + win];
				}
			}

			// ===========================================================
			// Perform mean substraction and normalization
			// ===========================================================             
			if (bDetrend){
				// Get mean
				for(win=0;win<nSamples;win++){
					avg1 += delayedLeft[win];
					avg2 += delayedRight[win];
				}
				avg1/=nSamples;
				avg2/=nSamples;

				if (bNorm){
					for(win=0;win<nSamples;win++) {
						// Auto-correlation (detrend)
						d1 +=(delayedLeft[win]-avg1)  * (delayedLeft[win]-avg1);
						d2 +=(delayedRight[win]-avg2) * (delayedRight[win]-avg2);

						// Cross-correlation (detrend)
						cc += (delayedLeft[win]-avg1) * (delayedRight[win]-avg2);
					}
					// Scale denominator
					d1=sqrt(d1);
					d2=sqrt(d2);
					// Normalization (added EPS to avoid NaN's if zero-padding is used to
					// extend the correlation range above the number of available samples)
					output[chan*nLags + delay + maxLag] = cc / (EPS + (double)(d1*d2));
				}
				else{
					for(win=0;win<nSamples;win++) {
						// Cross-correlation (detrend)
						cc += (delayedLeft[win]-avg1) * (delayedRight[win]-avg2);
					}
					output[chan*nLags + delay + maxLag] = cc;
				}
			}
			else{
				if (bNorm){
					for(win=0;win<nSamples;win++) {
						// Auto-correlation
						d1 +=(delayedLeft[win])  * (delayedLeft[win]);
						d2 +=(delayedRight[win]) * (delayedRight[win]);
						// Cross-correlation
						cc += (delayedLeft[win]) * (delayedRight[win]);
					}
					// Scale denominator
					d1=sqrt(d1);
					d2=sqrt(d2);
					// Normalization (added EPS to avoid NaN's if zero-padding is used to
					// extend the correlation range above the number of available samples)
					output[chan*nLags + delay + maxLag] = cc / (EPS + (double)(d1*d2));
				}
				else{
					for(win=0;win<nSamples;win++) {
						// Cross-correlation
						cc += (delayedLeft[win]) * (delayedRight[win]);
					}
					output[chan*nLags + delay + maxLag] = cc;
				}
			}
		}
	}
}   /* end mexFunction() */

