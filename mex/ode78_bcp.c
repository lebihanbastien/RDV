/*=========================================================================
 * ode78_bcp.c mex file to generate
 * a Runge-Kutta 7/8 integrator of the CR3BP vector field.
 * 
 * author:  BLB
 * year:    2015
 * version: 1.0
 *=======================================================================*/

//-------------------------------------------------------------------------
// Headers
//-------------------------------------------------------------------------
//Mex
#include "mex.h"

//Custom
#include "C/ode78.h"


//-------------------------------------------------------------------------
// The gateway function.
// The input must bet, in that order:
// 1. double t0, the initial time
// 2. double tf, the final time
// 3. double y0[6 or 42], the initial state
// 4. int    nvar = 6 or 42, the number of state variables
// 5. double mu, the cr3bp mass ratio
// 6. double omega0, the initial phase of the Sun
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //---------------------------------------------------------------------
    // Check for proper number of arguments
    //---------------------------------------------------------------------
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("custom:ode78_bcp:nrhs","6 inputs required.");
    }
    if(!(nlhs == 2 || nlhs == 4)) {
        mexErrMsgIdAndTxt("custom:ode78_bcp:nlhs",
       "2 or 4 outputs required: 2 if only the final state is desired, 4 if the state on a given time grid is also desired.");
    }
 
    //---------------------------------------------------------------------
    // Retrieve the variables:
    // 1. t0, the initial time
    // 2. tf, the final time
    // 3. y0, the initial state
    // 4. nvar, the number of state variables
    // 5. mu, the cr3bp mass ratio
    // 6. om0, the initial phase angle of the Sun
    //---------------------------------------------------------------------
    double t0  = mxGetScalar(prhs[0]);
    double tf  = mxGetScalar(prhs[1]);
    double *y0 = mxGetPr(prhs[2]);
    int nvar   = (int) mxGetScalar(prhs[3]);
    double mu  = mxGetScalar(prhs[4]);
    double om0 = mxGetScalar(prhs[4]);
    //---------------------------------------------------------------------
    // State will be stored on a given grid, if necessary
    //---------------------------------------------------------------------
    int nGrid = 1000;
    //Event states
    double **yv = (double **) calloc(nvar, sizeof(double*));
    for(int n = 0; n < nvar; n++) yv[n] = (double*) calloc(nGrid+1, sizeof(double));
    //Event time
    double *tv = (double*) calloc(nGrid+1, sizeof(double));
  
    //---------------------------------------------------------------------
    // BCP params
    //---------------------------------------------------------------------
    double *params = (double*) calloc(2, sizeof(double));
    params[0] = mu;
    params[1] = om0;
    
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    double t, y[nvar];
    if(nlhs == 2) ode78_bcp(&t, y, y0, t0, tf, nvar, params);
    if(nlhs == 4) ode78_bcp_vec(&t, y, tv, yv, nGrid, y0, t0, tf, nvar, params);
        
    //---------------------------------------------------------------------
    // Output: the final state
    //---------------------------------------------------------------------
    //Create the output matrices
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);           //t
    plhs[1] = mxCreateDoubleMatrix(1, nvar, mxREAL);        //y
    

    //Get a pointer to the real data in the output
    double *tout  = mxGetPr(plhs[0]);
    double *yout  = mxGetPr(plhs[1]);
    
    //Store the final state
    for(int i = 0; i < nvar; i++) yout[i] = y[i];
    tout[0] = t;
    
    //---------------------------------------------------------------------
    // Output: the state on a time grid, if necessary
    //---------------------------------------------------------------------
    if(nlhs == 4)
    {
        //Create the output matrices
        plhs[2] = mxCreateDoubleMatrix(nGrid+1, 1, mxREAL);     //tv
        plhs[3] = mxCreateDoubleMatrix(nGrid+1, nvar, mxREAL);  //yv
           
        //Get a pointer to the real data in the output
        double *tvout = mxGetPr(plhs[2]);
        double *yvout = mxGetPr(plhs[3]);
        
        //Store the state on the grid [0,..., nGrid]
        int indix = 0;
        for(int i = 0; i < nvar; i++)
        {
            for(int k = 0; k <= nGrid; k++)
            {
                yvout[indix++] = yv[i][k];
            }
        }
        for(int k = 0; k <= nGrid; k++) tvout[k] = tv[k];

    }
    
    //---------------------------------------------------------------------
    // Free memory
    //---------------------------------------------------------------------
    free(tv);
    free(params);
    for (int i=0; i<=nvar; i++) free(yv[i]); free(yv);
}