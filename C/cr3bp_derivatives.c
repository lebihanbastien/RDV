#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

//Custom
#include "cr3bp_derivatives.h"
#include "custom.h"

/**
* \brief compute the derivatives of  the state variables [x, y, z, xp, yp, zp]
* \param t: the current time (unused since the system is autonomous)
* \param y: the current state (dim 6)
* \param f: the derivatives to update
* \param: *params: a pointer towards the integration paramater mu
**/
int cr3bp_derivatives_6 (double t, const double y[], double f[], void *params)
{
    double mu = *(double *)params;

    //-------------------------------------------------------------------------------
    //Update first & second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //-------------------------------------------------------------------------------
    double dU[3];

    dU[0] = (mu*(2*mu + 2*y[0] - 2))/(2*(pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0))) - ((2*mu + 2*y[0])*(mu - 1))/(2*pow(pow(mu + y[0],2.0) + pow(y[1],2.0) + pow(y[2],2.0),3.0/2.0)) - y[0];
    dU[1] = (mu*y[1])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - y[1] - (y[1]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),3/2.0);
    dU[2] = (mu*y[2])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (y[2]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),3/2.0);

    //-------------------------------------------------------------------------------
    //Phase space derivatives
    //-------------------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = -dU[0] + 2*y[4];
    f[4] = -dU[1] - 2*y[3];
    f[5] = -dU[2];


    return GSL_SUCCESS;
}

/**
* \brief compute the derivatives of both the state variables [x, y, z, xp, yp, zp] and the STM in a single output f of 42 dim
* \param t: the current time (unused since the system is autonomous)
* \param y: the current state (dim 42)
* \param f: the derivatives to update
* \param: *params: a pointer towards the integration paramater mu
**/
int cr3bp_derivatives_42 (double t, const double y[], double f[], void *params)
{

    double mu = *(double *)params;
    int i,j;

    //-------------------------------------------------------------------------------
    //Update first & second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //-------------------------------------------------------------------------------
    double dU[3];
    double dU2[3][3];

    dU[0] = (mu*(2*mu + 2*y[0] - 2))/(2*(pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0))) - ((2*mu + 2*y[0])*(mu - 1))/(2*pow(pow(mu + y[0],2.0) + pow(y[1],2.0) + pow(y[2],2.0),3.0/2.0)) - y[0];
    dU[1] = (mu*y[1])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - y[1] - (y[1]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),3/2.0);
    dU[2] = (mu*y[2])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (y[2]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),3/2.0);

    dU2[0][0] =  mu/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (mu - 1)/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2), 3/2.0) - (3*mu*pow(2*mu + 2*y[0] - 2,2))/(4*pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0)) + (3*pow(2*mu + 2*y[0],2)*(mu - 1))/(4*pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0)) - 1;
    dU2[0][1] = (3*y[1]*(2*mu + 2*y[0])*(mu - 1))/(2*pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0)) - (3*mu*y[1]*(2*mu + 2*y[0] - 2))/(2*pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0));
    dU2[0][2] = (3*y[2]*(2*mu + 2*y[0])*(mu - 1))/(2*pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0)) - (3*mu*y[2]*(2*mu + 2*y[0] - 2))/(2*pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0));
    dU2[1][0] = dU2[0][1];
    dU2[1][1] = mu/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (mu - 1)/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2), 3/2.0) + (3*pow(y[1],2)*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0) - (3*mu*pow(y[1],2))/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0) - 1;
    dU2[1][2] = (3*y[1]*y[2]*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0) - (3*mu*y[1]*y[2])/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0);
    dU2[2][0] = dU2[0][2];
    dU2[2][1] = dU2[1][2];
    dU2[2][2] = mu/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),3/2.0) - (mu - 1)/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2), 3/2.0) + (3*pow(y[2],2)*(mu - 1))/pow(pow(mu + y[0],2) + pow(y[1],2) + pow(y[2],2),5/2.0) - (3*mu*pow(y[2],2))/pow(pow(mu + y[0] - 1,2) + pow(y[1],2) + pow(y[2],2),5/2.0);

    //GSL version
    gsl_matrix *dU2_gsl = gsl_matrix_alloc(3,3);
    for(i=0; i<3 ; i++)
        for(j=0; j<3; j++) gsl_matrix_set(dU2_gsl, i, j, dU2[i][j]);
    gsl_matrix_scale (dU2_gsl, -1); //take the opposite of dU2 for later concatenation

    //-------------------------------------------------------------------------------
    //Phase space derivatives
    //-------------------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = -dU[0] + 2*y[4];
    f[4] = -dU[1] - 2*y[3];
    f[5] = -dU[2];


    //-------------------------------------------------------------------------------
    //STM derivative
    //-------------------------------------------------------------------------------
    //STM is updated
    gsl_matrix *STM = gsl_matrix_alloc(6,6);
    custom_vectorToMatrix(STM, (double *)y, 6, 6, 6);

    //Trivial matrices
    gsl_matrix *eye3_gsl = gsl_matrix_alloc(3,3);
    gsl_matrix *z3_gsl = gsl_matrix_alloc(3,3);
    gsl_matrix *Omega2_gsl = gsl_matrix_alloc(3,3);

    gsl_matrix_set_zero (z3_gsl);
    gsl_matrix_set_identity (eye3_gsl);
    gsl_matrix_set_zero(Omega2_gsl);

    gsl_matrix_set(Omega2_gsl, 0, 1, 2);
    gsl_matrix_set(Omega2_gsl, 1, 0, -2);

    // Building Df
    //----------------------------------------------------------------------------
    //
    //Matrix concatenation:
    //
    //      |  O          I3  |
    //Df =  |  -dU2    2Omega |
    gsl_matrix *Df_gsl = gsl_matrix_alloc(6,6);

    gsl_matrix_view Df_ul = gsl_matrix_submatrix (Df_gsl , 0 , 0 , 3 , 3 );
    gsl_matrix_view Df_ur = gsl_matrix_submatrix (Df_gsl , 0 , 3 , 3 , 3 );
    gsl_matrix_view Df_ll = gsl_matrix_submatrix (Df_gsl , 3 , 0 , 3 , 3 );
    gsl_matrix_view Df_lr = gsl_matrix_submatrix (Df_gsl , 3 , 3 , 3 , 3 );

    gsl_matrix_memcpy( &Df_ul.matrix, z3_gsl);
    gsl_matrix_memcpy( &Df_ur.matrix, eye3_gsl);
    gsl_matrix_memcpy( &Df_ll.matrix, dU2_gsl);
    gsl_matrix_memcpy( &Df_lr.matrix, Omega2_gsl);

    //Matrix product
    //----------------------------------------------------------------------------
    gsl_matrix *STMp = gsl_matrix_alloc(6,6);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Df_gsl, STM,0.0, STMp);

    //----------------------------------------------------------------------------
    //Update f (from row 6 to 41)
    //----------------------------------------------------------------------------
    custom_matrixToVector(f, STMp, 6, 6, 6);


    //Memory release
    //----------------------------------------------------------------------------
    gsl_matrix_free(dU2_gsl);
    gsl_matrix_free(STM);
    gsl_matrix_free(eye3_gsl);
    gsl_matrix_free(z3_gsl);
    gsl_matrix_free(Omega2_gsl);
    gsl_matrix_free(STMp);
    gsl_matrix_free(Df_gsl);

    return GSL_SUCCESS;


}



/**
* \brief compute the derivatives of  the state variables [x, y, z, xp, yp, zp], for the bicircular problem
* \param t: the current time (unused since the system is autonomous)
* \param y: the current state (dim 6)
* \param f: the derivatives to update
* \param: *params: a pointer towards the integration paramater mu
**/
int bcp_derivatives_6 (double t, const double y[], double f[], void *params)
{
    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    double *p = (double *)params;
    double mu     = p[0];    //crtbp mass ratio
    double omega0 = p[1];    //initial phase of the fourth body
    
    //printf("mu = %5.5f, omega0 = %5.5f\n", mu, omega0);
    
    //Sun parameters
    double ms = 328900.54; 
    double as = 388.81114;
    double omegaS = -0.925195985520347;

    //Current Sun phase angle
    double theta = omega0 + omegaS*t; 
    
    //-------------------------------------------------------------------------------
    //Update first & second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //-------------------------------------------------------------------------------
    double r1 = sqrt((y[0]+mu)*(y[0]+mu) + y[1]*y[1] + y[2]*y[2]);
    double r2 = sqrt((y[0]-1+mu)*(y[0]-1+mu) + y[1]*y[1] + y[2]*y[2]);
    double rs = sqrt((y[0]-as*cos(theta))*(y[0]-as*cos(theta)) + (y[1]-as*sin(theta))*(y[1]-as*sin(theta)) + y[2]*y[2]);
    
    double dU[3];
    dU[0] = y[0] - (1.0-mu)/pow(r1, 3.0)*(y[0]+mu) - mu/pow(r2, 3.0)*(y[0]-1+mu) - ms/pow(rs, 3.0)*(y[0]-as*cos(theta)) - ms/pow(as, 2.0)*cos(theta);
    dU[1] = y[1] - (1.0-mu)/pow(r1, 3.0)*y[1]      - mu/pow(r2, 3.0)*y[1]        - ms/pow(rs, 3.0)*(y[1]-as*sin(theta)) - ms/pow(as, 2.0)*sin(theta);
    dU[2] =      - (1.0-mu)/pow(r1, 3.0)*y[2]      - mu/pow(r2, 3.0)*y[2]        - ms/pow(rs, 3.0)*y[2];
  
    //-------------------------------------------------------------------------------
    //Phase space derivatives: careful, change of sign wrt CR3BP case!!
    //-------------------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = +dU[0] + 2*y[4];
    f[4] = +dU[1] - 2*y[3];
    f[5] = +dU[2];


    return GSL_SUCCESS;
}



/**
* \brief compute the derivatives of both the state variables [x, y, z, xp, yp, zp] and the STM in a single output f of 42 dim, for the bicircular problem
* \param t: the current time (unused since the system is autonomous)
* \param y: the current state (dim 42)
* \param f: the derivatives to update
* \param: *params: a pointer towards the integration paramater mu
**/
int bcp_derivatives_42 (double t, const double y[], double f[], void *params)
{
    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    double *p = (double *)params;
    double mu     = p[0];    //crtbp mass ratio
    double omega0 = p[1];    //initial phase of the fourth body
    
    //Sun parameters
    double ms = 328900.54; 
    double as = 388.81114;
    double omegaS = -0.925195985520347;

    //Current Sun phase angle
    double theta = omega0 + omegaS*t; 
    
    int i,j;

    //-------------------------------------------------------------------------------
    //Update first derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //-------------------------------------------------------------------------------
    double r1 = sqrt((y[0]+mu)*(y[0]+mu) + y[1]*y[1] + y[2]*y[2]);
    double r2 = sqrt((y[0]-1+mu)*(y[0]-1+mu) + y[1]*y[1] + y[2]*y[2]);
    double rs = sqrt((y[0]-as*cos(theta))*(y[0]-as*cos(theta)) + (y[1]-as*sin(theta))*(y[1]-as*sin(theta)) + y[2]*y[2]);
      
    double dU[3];
    dU[0] = y[0] - (1.0-mu)/pow(r1, 3.0)*(y[0]+mu) - mu/pow(r2, 3.0)*(y[0]-1+mu) - ms/pow(rs, 3.0)*(y[0]-as*cos(theta)) - ms/pow(as, 2.0)*cos(theta);
    dU[1] = y[1] - (1.0-mu)/pow(r1, 3.0)*y[1]      - mu/pow(r2, 3.0)*y[1]        - ms/pow(rs, 3.0)*(y[1]-as*sin(theta)) - ms/pow(as, 2.0)*sin(theta);
    dU[2] =      - (1.0-mu)/pow(r1, 3.0)*y[2]      - mu/pow(r2, 3.0)*y[2]        - ms/pow(rs, 3.0)*y[2];

    //-------------------------------------------------------------------------------
    //Update second derivatives of the potential \bar{U} (cf Koon et al. 2006)
    //-------------------------------------------------------------------------------
    double dU2[3][3];
    //dUx/dx,y,z
    dU2[0][0] = 1 + 3.0*(1.0-mu)/pow(r1, 5.0)*(y[0]+mu)*(y[0]+mu)                 - (1.0-mu)/pow(r1, 3.0)
                  + 3.0*mu/pow(r2, 5.0)*(y[0]-1.0+mu)*(y[0]-1.0+mu)               - mu/pow(r2, 3.0)
                  + 3.0*ms/pow(rs, 5.0)*(y[0]-as*cos(theta))*(y[0]-as*cos(theta)) - ms/pow(rs, 3.0);
    
    dU2[0][1] = 0 + 3.0*(1.0-mu)/pow(r1, 5.0)*(y[0]+mu)*y[1]
                  + 3.0*mu/pow(r2, 5.0)*(y[0]-1.0+mu)*y[1]
                  + 3.0*ms/pow(rs, 5.0)*(y[0]-as*cos(theta))*(y[1]-as*sin(theta));
    
    dU2[0][2] = 0 + 3.0*(1.0-mu)/pow(r1, 5.0)*(y[0]+mu)*y[2]
                  + 3.0*mu/pow(r2, 5.0)*(y[0]-1.0+mu)*y[2]
                  + 3.0*ms/pow(rs, 5.0)*(y[0]-as*cos(theta))*y[2];
    
    //dUy/dx,y,z
    dU2[1][0] = dU2[0][1];
    
    dU2[1][1] = 1 + 3.0*(1.0-mu)/pow(r1, 5.0)*y[1]*y[1]                           - (1.0-mu)/pow(r1, 3.0)
                  + 3.0*mu/pow(r2, 5.0)*y[1]*y[1]                                 - mu/pow(r2, 3.0)
                  + 3.0*ms/pow(rs, 5.0)*(y[1]-as*sin(theta))*(y[1]-as*sin(theta)) - ms/pow(rs, 3.0);
    
    dU2[1][2] = 0 + 3.0*(1.0-mu)/pow(r1, 5.0)*y[1]*y[2]
                  + 3.0*mu/pow(r2, 5.0)*y[1]*y[2]
                  + 3.0*ms/pow(rs, 5.0)*(y[1]-as*sin(theta))*y[2];
   
   //dUz/dx,y,z             
   dU2[2][0]  = dU2[0][2];
   dU2[2][1]  = dU2[1][2];
   dU2[2][2]  = 0 + 3.0*(1.0-mu)/pow(r1, 5.0)*y[2]*y[2] - (1.0-mu)/pow(r1, 3.0)
                  + 3.0*mu/pow(r2, 5.0)*y[2]*y[2]       - mu/pow(r2, 3.0)
                  + 3.0*ms/pow(rs, 5.0)*y[2]*y[2]       - ms/pow(rs, 3.0);    
    
    //GSL version
    gsl_matrix *dU2_gsl = gsl_matrix_alloc(3,3);
    for(i=0; i<3 ; i++)
        for(j=0; j<3; j++) gsl_matrix_set(dU2_gsl, i, j, dU2[i][j]);
    //gsl_matrix_scale (dU2_gsl, -1); //do NOT take the opposite of dU2 here!!
    

    //-------------------------------------------------------------------------------
    //Phase space derivatives: careful, change of sign wrt CR3BP case!!
    //-------------------------------------------------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = dU[0] + 2*y[4];
    f[4] = dU[1] - 2*y[3];
    f[5] = dU[2];


    //-------------------------------------------------------------------------------
    //STM derivative
    //-------------------------------------------------------------------------------
    //STM is updated
    gsl_matrix *STM = gsl_matrix_alloc(6,6);
    custom_vectorToMatrix(STM, (double *)y, 6, 6, 6);

    //Trivial matrices
    gsl_matrix *eye3_gsl = gsl_matrix_alloc(3,3);
    gsl_matrix *z3_gsl = gsl_matrix_alloc(3,3);
    gsl_matrix *Omega2_gsl = gsl_matrix_alloc(3,3);

    gsl_matrix_set_zero (z3_gsl);
    gsl_matrix_set_identity (eye3_gsl);
    gsl_matrix_set_zero(Omega2_gsl);

    gsl_matrix_set(Omega2_gsl, 0, 1, 2);
    gsl_matrix_set(Omega2_gsl, 1, 0, -2);

    // Building Df
    //----------------------------------------------------------------------------
    //
    //Matrix concatenation:
    //
    //      |  O          I3  |
    //Df =  |  -dU2    2Omega |
    gsl_matrix *Df_gsl = gsl_matrix_alloc(6,6);

    gsl_matrix_view Df_ul = gsl_matrix_submatrix (Df_gsl , 0 , 0 , 3 , 3 );
    gsl_matrix_view Df_ur = gsl_matrix_submatrix (Df_gsl , 0 , 3 , 3 , 3 );
    gsl_matrix_view Df_ll = gsl_matrix_submatrix (Df_gsl , 3 , 0 , 3 , 3 );
    gsl_matrix_view Df_lr = gsl_matrix_submatrix (Df_gsl , 3 , 3 , 3 , 3 );

    gsl_matrix_memcpy( &Df_ul.matrix, z3_gsl);
    gsl_matrix_memcpy( &Df_ur.matrix, eye3_gsl);
    gsl_matrix_memcpy( &Df_ll.matrix, dU2_gsl);
    gsl_matrix_memcpy( &Df_lr.matrix, Omega2_gsl);

    //Matrix product
    //----------------------------------------------------------------------------
    gsl_matrix *STMp = gsl_matrix_alloc(6,6);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Df_gsl, STM,0.0, STMp);

    //----------------------------------------------------------------------------
    //Update f (from row 6 to 41)
    //----------------------------------------------------------------------------
    custom_matrixToVector(f, STMp, 6, 6, 6);


    //Memory release
    //----------------------------------------------------------------------------
    gsl_matrix_free(dU2_gsl);
    gsl_matrix_free(STM);
    gsl_matrix_free(eye3_gsl);
    gsl_matrix_free(z3_gsl);
    gsl_matrix_free(Omega2_gsl);
    gsl_matrix_free(STMp);
    gsl_matrix_free(Df_gsl);

    return GSL_SUCCESS;


}
