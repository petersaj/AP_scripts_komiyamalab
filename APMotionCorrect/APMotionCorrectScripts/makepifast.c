/*=================================================================
 *
 * YPRIME.C	Sample .MEX file corresponding to YPRIME.M
 *	        Solves simple 3 body orbit problem 
 *
 * The calling syntax is:
 *
 *		[yp] = yprime(t, y)
 *
 *  You may also want to look at the corresponding M-code, yprime.m.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2006 The MathWorks, Inc.
 *
 *=================================================================*/
/* $Revision: 1.10.6.4 $ */
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	MAT_IN	prhs[0]

/* Output Arguments */

#define	MEAN_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    mxArray *xData,*x2Data;
    double *xVal,*x2Val;
    double *outArray;
    double *y,**x;  
    int rowLen,colLen,i,j;
    double avg=0;
    int F,G;
    double mu;
    
    /* Check for proper number of arguments */
   
    xData=prhs[0];
    x2Data=prhs[1];

    //Get matrix x
    xVal = mxGetPr(xData);
    x2Val = mxGetPr(x2Data);
    rowLen = mxGetN(xData);
    colLen = mxGetM(xData);
    plhs[0] = mxCreateDoubleMatrix(rowLen, colLen, mxREAL);
    outArray = mxGetPr(plhs[0]);
    
    
    for(i=0;i<rowLen;i++)
    {
        for(j=0;j<colLen;j++)
        {
            //avg += xVal[(i*colLen)+j];
   
            //outArray[(i*colLen)+j]=xVal[(i*colLen)+j]*x2Val[j];
            outArray[(i*colLen)+j]=xVal[(i*colLen)+j]+x2Val[j];
            
        }
      
       
    }

    return;

}


