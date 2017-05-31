//Custom
#include "define_env.h"
#include "cr3bp_derivatives.h"
#include "custom_odezero.h"
#include "custom_ode.h"
#include "differential_correction.h"
#include "user_parameters.h"

int differential_correction(CR3BP cr3bp,
                            double ystart[],
                            double yhalf[],
                            double *tcross,
                            double t1,
                            custom_ode_structure *ode_s,
                            diff_corr_parameters diff_corr)
{
    //Current error on -xf and -zf
    double dyf[2];
    double y[42];
    int i,j;
    int iter = 0;
    int status;
    do
    {
        //Update the starting point
        for(i=0; i<42; i++) y[i] = ystart[i];

        //Integration until y=0 is reached
        status = custom_odezero(y, tcross, t1, ode_s);

        if(status != GSL_SUCCESS){
            printf("WARNING: odezero failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //Af matrix
        double Af[2][2];
        switch(diff_corr.diff_corr_type)
        {
        case Z0_FIXED:
            Af[0][0]= y[5+19];  //Phi41
            Af[0][1]= y[5+23];  //Phi45
            Af[1][0]= y[5+31];  //Phi61
            Af[1][1]= y[5+35];  //Phi65
            break;

        case X0_FIXED:
            Af[0][0]= y[5+21];  //Phi43
            Af[0][1]= y[5+23];  //Phi45
            Af[1][0]= y[5+33];  //Phi63
            Af[1][1]= y[5+35];  //Phi65
            break;

            default: //default is Z0_FIXED

            Af[0][0]= y[5+19];  //Phi41
            Af[0][1]= y[5+23];  //Phi45
            Af[1][0]= y[5+31];  //Phi61
            Af[1][1]= y[5+35];  //Phi65
            break;
        }

        //Bf vector
        double Bf[2];
        switch(diff_corr.diff_corr_type)
        {
        case Z0_FIXED:
            Bf[0] = y[5+7];  //Phi21
            Bf[1] = y[5+11]; //Phi25
            break;

        case X0_FIXED:
            Bf[0] = y[5+9];  //Phi23
            Bf[1] = y[5+11]; //Phi25
            break;

            default: //default is Z0_FIXED
            Bf[0] = y[5+7];  //Phi21
            Bf[1] = y[5+11]; //Phi25
            break;
        }


        //State subvector
        double y_state[6];
        for(i =0; i<6; i++) y_state[i] = y[i];

        //Derivative of y_state @ y=0 in y_state_p (tcross is not used)
        double y_state_p[6];
        cr3bp_derivatives_6 (*tcross, y_state, y_state_p, &cr3bp.mu);

        //ppf vector
        double ppf[2];
        ppf[0] = y_state_p[3]/y_state[4];
        ppf[1] = y_state_p[5]/y_state[4];

        //Af_new
        double Af_new[2][2];
        for(i=0; i<2; i++)
        {
            for(j=0; j<2; j++) Af_new[i][j] = Af[i][j] - ppf[i]*Bf[j];
        }

        //Inversion of the error
        dyf[0] = -y_state[3];
        dyf[1] = -y_state[5];
        double dy0[2];
        dy0[1] = (Af_new[0][0]*dyf[1] - Af_new[1][0]*dyf[0])/(Af_new[1][1]*Af_new[0][0]-Af_new[0][1]*Af_new[1][0]);
        dy0[0] = (dyf[0] - Af_new[0][1]*dy0[1])/Af_new[0][0];

        //Updating the state
        switch(diff_corr.diff_corr_type)
        {
        case Z0_FIXED:
            ystart[0] = ystart[0] + dy0[0];  //x0
            ystart[4] = ystart[4] + dy0[1];  //dy0
            break;

        case X0_FIXED:
            ystart[2] = ystart[2] + dy0[0];  //z0
            ystart[4] = ystart[4] + dy0[1];  //dy0
            break;
        }

    }
    while(fabs(dyf[0])> diff_corr.eps_diff && fabs(dyf[1])> diff_corr.eps_diff && (++iter) < 50);

    if(iter>=50){
        printf("WARNING: number of iter max exceeded during differential correction. Premature ending\n");
        return GSL_FAILURE;
    }

    //Update the half-period state
    for(i=0;i<42;i++) yhalf[i] = y[i];

    return GSL_SUCCESS;
}
