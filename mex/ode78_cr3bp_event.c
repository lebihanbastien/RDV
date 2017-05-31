/*=========================================================================
 * ode78_cr3bp.c mex file to generate
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
#include "C/custom_ode.h"
#include "C/custom_odezero.h"
#include "C/cr3bp_derivatives.h"

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

//-------------------------------------------------------------------------
// C function
//-------------------------------------------------------------------------
void ode78_cr3bp_event(double *t,                    //current time
                       double *y,                    //current state
                       double **ye,                  //event states
                       double *te,                   //event times
                       double const *y0,             //initial condition
                       double const t0,              //initial time
                       double tf,                    //final time
                       int nvar,                     //number of state variables
                       struct value_params *val_par, //event structure
                       int *nEvents,                 //the effective number of events found
                       double *mu)                   //cr3bp mass ratio
{
    //---------------------------------------------------------------------
    // Initialize the integration structures
    //---------------------------------------------------------------------
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    
    
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method

    //General ode structure
    custom_ode_structure ode_s;
    switch(nvar)
    {
        case 6:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 6, 1e-6, cr3bp_derivatives_6, NULL, mu);
            break;
        }
        case 42:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 42, 1e-6, cr3bp_derivatives_42, NULL, mu);
            break;
        }
        default:   //if nvar !=6 && nvar != 42, return without integration
        {
            mexPrintf("ode78_cr3bp. Error: wrong number nvar of state variables. nvar must 6 or 42. return.");
            return;
        }
    
    }
    
    
    //---------------------------------------------------------------------
    // Initalization of the state
    //---------------------------------------------------------------------
    *t = t0;
    for(int i = 0; i < nvar; i++) y[i] = y0[i];   
    
    //---------------------------------------------------------------------
    // Integration until t = t1, and saving the events
    //---------------------------------------------------------------------
    struct value_function fvalue;
    fvalue.val_par = val_par;
    
    switch(val_par->type)
    {
        case 'X':
        case 'Y':
        case 'Z':
        {
            fvalue.value   = &linear_intersection; 
            break;
        }
        
        case 'A':
        {
            fvalue.value   = &angle_intersection; 
            break;
        }
        
        case 'F':
        {
            fvalue.value   = &null_flight_path_angle; 
            break;
        }
        
        default:
        {
            mexPrintf("Bad event type. return.");
            return;
        }
    }
    
    
    //Starting in the right direction
    ode_s.d->h = (tf>t0) ? fabs(ode_s.d->h) : -fabs(ode_s.d->h);
    ode_s.h = (tf>t0) ? fabs(ode_s.h) : -fabs(ode_s.h);
    
    //Apply custom_odezero_2
    *nEvents = custom_odezero_2(y, ye, te, tf, &ode_s, fvalue);
    
    // custom_odezero_2 return the last indix rather than the number of events
    // So we need to add one event
    *nEvents = *nEvents+1;
}


//-------------------------------------------------------------------------
// The gateway function.
// The inputs must be, in that order:
// 1. double t0, the initial time
// 2. double tf, the final time
// 3. double y0[6 or 42], the initial state
// 4. int    nvar = 6 or 42, the number of state variables
// 5. double mu, the cr3bp mass ratio
// 6. struct struct_event the event structure 
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //---------------------------------------------------------------------
    // Check for proper number of arguments
    //---------------------------------------------------------------------
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("custom:ode78_cr3bp_events:nrhs","6 inputs required.");
    }
    if(!(nlhs == 2 || nlhs == 4)) {
        mexErrMsgIdAndTxt("custom:ode78_cr3bp_events:nlhs",
       "2 or 4 outputs required: 2 if only the final state is desired, 4 if the state on a given time grid is also desired.");
    }
 
    //---------------------------------------------------------------------
    // Retrieve the variables:
    // 1. t0, the initial time
    // 2. tf, the final time
    // 3. y0, the initial state
    // 4. nvar, the number of state variables
    // 5. mu, the cr3bp mass ratio
    // 6. struct struct_event the event structure 
    //---------------------------------------------------------------------
    double t0  = mxGetScalar(prhs[0]);
    double tf  = mxGetScalar(prhs[1]);
    double *y0 = mxGetPr(prhs[2]);
    int nvar   = (int) mxGetScalar(prhs[3]);
    double mu  = mxGetScalar(prhs[4]);
    
    //---------------------------------------------------------------------
    // particular case of the event structure:
    //---------------------------------------------------------------------
    // Get input arguments
    mwSize nfields = mxGetNumberOfFields(prhs[5]);           //number of fields
    mwSize NStructElems = mxGetNumberOfElements(prhs[5]);    //number of elements (must be 1)
    
    if(NStructElems != 1)
    {
        mexErrMsgIdAndTxt("custom:ode78_cr3bp_events:NStructElems",
                "Only one element is required in the event structure.");
    }
    
    
    if(nfields != 7)
    {
        printf(" nfields = %d", nfields);
        mexErrMsgIdAndTxt("custom:ode78_cr3bp_events:nfields",
                "The event structure must contain 7 fields: type, isterminal, max_events, direction, dim, value, and center.");
    }
    
    //---------------------------------------------------------------------
    //Get the different fields in the event structure val_par = prhs[5]:
    // val_par.max_events
    // val_par.direction
    // val_par.dim
    // val_par.value
    //---------------------------------------------------------------------
    mxArray *tmp;
    int max_events, direction, dim;
    char *type;
    double value;
    double *center = (double*) calloc(3, sizeof(double));

    //0. Type
    tmp = mxGetField(prhs[5], 0, "type");
    type = (char*) mxGetChars(tmp);
    //mexPrintf("type = %s\n", type);
    
    //1. val_par.max_events
    tmp = mxGetField(prhs[5], 0, "max_events");
    max_events = (int) mxGetScalar(tmp);
    //mexPrintf("max_events = %d\n", max_events); 
    
    //2. val_par.direction
    tmp = mxGetField(prhs[5], 0, "direction");
    direction = (int) mxGetScalar(tmp);
    //mexPrintf("direction = %d\n", direction); 
    
    //3. val_par.dim
    tmp = mxGetField(prhs[5], 0, "dim");
    dim = (int) mxGetScalar(tmp)-1;
    //mexPrintf("dim = %d\n", dim); 
    
    //4. val_par.value
    tmp = mxGetField(prhs[5], 0, "value");
    value = mxGetScalar(tmp);
    //mexPrintf("value = %5.5f\n", value); 
    
    //5. val_par.center
    tmp = mxGetField(prhs[5], 0, "center");
    center = mxGetPr(tmp);
    //mexPrintf("center = (%5.5f, %5.5f, %5.5f)\n", center[0], center[1], center[2]);
    
    
    
    //---------------------------------------------------------------------
    // Create the corresponding C structure
    //---------------------------------------------------------------------
    struct value_params val_par;
    val_par.max_events = max_events;
    val_par.direction  = direction;
    val_par.dim        = dim;
    val_par.value      = value;
    val_par.center     = center;
    val_par.type       = (int) *type;
    
    //---------------------------------------------------------------------
    // Number of events
    //---------------------------------------------------------------------
    int MAX_EVENTS = 50;  //maximum allowed is 50
    int nEv; //actual number of events stored at the end of the computation
    
    if(max_events > MAX_EVENTS)
    {
        mexErrMsgIdAndTxt("custom:ode78_cr3bp_events:maxEv","Number of required events must be < 50.");
    }
    
    
    //---------------------------------------------------------------------
    // Initialization
    //---------------------------------------------------------------------
    //Current state & time
    double t, y[nvar];
    //Event states
    double **ye = (double **) calloc(nvar, sizeof(double*));
    for(int n = 0; n < nvar; n++) ye[n] = (double*) calloc(MAX_EVENTS+1, sizeof(double));
    //Event time
    double *te = (double*) calloc(MAX_EVENTS+1, sizeof(double));
    
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78_cr3bp_event(&t, y, ye, te, y0, t0, tf, nvar, &val_par, &nEv, &mu);
    
    //print the number of events found
    //mexPrintf("ode78_cr3bp_event. %d were found.\n", nEv);
    
    //---------------------------------------------------------------------
    // If it is desired, store the state on a time grid up to te[end]
    //---------------------------------------------------------------------
    int nGrid = 1000;
    //Event states
    double **yv = (double **) calloc(nvar, sizeof(double*));
    for(int n = 0; n < nvar; n++) yv[n] = (double*) calloc(nGrid+1, sizeof(double));
    //Event time
    double *tv = (double*) calloc(nGrid+1, sizeof(double));
    //New final time
    tf = te[nEv-1];
    //Integration, if necessary
    if(nlhs == 4)
    {
        ode78_cr3bp_vec(&t, y, tv, yv, nGrid, y0, t0, tf, nvar, &mu);
    }
    //---------------------------------------------------------------------
    // Output
    //---------------------------------------------------------------------
    // create the output matrices
    plhs[0] = mxCreateDoubleMatrix(nEv, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nEv, nvar, mxREAL);
    
    //get a pointer to the real data in the output matrix
    double *tout = mxGetPr(plhs[0]);
    double *yout = mxGetPr(plhs[1]);
    
    
    //Store the event states
    int indix = 0;
    for(int i = 0; i < nvar; i++)
    {
        for(int k = 0; k < nEv; k++) 
        {
            yout[indix++] = ye[i][k];  
        }       
    }
    for(int k = 0; k < nEv; k++) tout[k] = te[k]; 
    
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

        //Print the gap between the last event and the grid computation
        //         mexPrintf("ode78_cr3bp_event. Discrepancy between last event and the grid computation\n");
        //         mexPrintf("ye    yvec    gap\n");
        //         for(int i = 0; i < nvar; i++)
        //         {
        //             mexPrintf("%5.5f     %5.5f      %5.5e \n", ye[i][nEv-1], yv[i][nGrid], ye[i][nEv-1]-yv[i][nGrid]);
        //         }
    }
    
    //---------------------------------------------------------------------
    // Free memory
    //---------------------------------------------------------------------
    free(tv);
    free(te);
    for (int i=0; i<=nvar; i++) free(yv[i]); free(yv);
    for (int i=0; i<=nvar; i++) free(ye[i]); free(ye);
}