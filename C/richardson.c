#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//Custom
#include "richardson.h"
#include "custom.h"
#include "orbits.h"

void third_order_halo_orbit(Orbit halo, double t, double yv_li[])
{
    int i;
    //---------------------------------------------------------------------------------------------------------------------------
    //Init
    //---------------------------------------------------------------------------------------------------------------------------
    double gamma_i = halo.Li.gamma_i;
    double mu = halo.cr3bp.mu;
    double dm = halo.dm;

    //Richardson coefficients
    RichardsonCoeffs RC;
    richardson_coefficients(halo.cr3bp, halo.Li.number, &RC);

    //X and Z-amplitude in Li-frame
    double Az_li = halo.Az_R/gamma_i;
    double Ax_li = sqrt(-(RC.l2*pow(Az_li,2.0) + RC.Delta)/RC.l1);

    //Phases
    double Phi = 0.0; //arbitrary
    //double Psi = m*M_PI/2;

    //nu parameter
    double nu = 1+RC.s1*pow(Ax_li,2.0)+RC.s2*pow(Az_li,2.0);

    //Third order estimate of the orbital period
    //double Tto = (2*M_PI)/(RC.omega_p*nu);

    //State estimate @ t  in the Li-frame
    double tau1 = RC.omega_p*nu*t + Phi;
    double dtau1dt = RC.omega_p*nu;

    yv_li[0] = RC.a21*pow(Ax_li,2.0) + RC.a22*pow(Az_li,2.0) - Ax_li*cos(tau1) + (RC.a23*pow(Ax_li,2.0) - RC.a24*pow(Az_li,2.0))*cos(2*tau1) + (RC.a31*pow(Ax_li,3.0) - RC.a32*Ax_li*pow(Az_li,2.0))*cos(3*tau1);
    yv_li[1] = RC.kappa*Ax_li*sin(tau1) + (RC.b21*pow(Ax_li,2.0) - RC.b22*pow(Az_li,2.0))*sin(2*tau1) + (RC.b31*pow(Ax_li,3.0) - RC.b32*Ax_li*pow(Az_li,2.0))*sin(3*tau1);
    yv_li[2] = dm*Az_li*cos(tau1) + dm*RC.d21*Ax_li*Az_li*(cos(2*tau1)-3) + dm*(RC.d32*Az_li*pow(Ax_li,2.0) - RC.d31*pow(Az_li,3.0))*cos(3*tau1);

    yv_li[3] = dtau1dt*( Ax_li*sin(tau1) - 2*(RC.a23*pow(Ax_li,2.0) - RC.a24*pow(Az_li,2.0))*sin(2*tau1) -3*(RC.a31*pow(Ax_li,3.0) - RC.a32*Ax_li*pow(Az_li,2.0))*sin(3*tau1)  );
    yv_li[4] = dtau1dt*( RC.kappa*Ax_li*cos(tau1) + 2*(RC.b21*pow(Ax_li,2.0) - RC.b22*pow(Az_li,2.0))*cos(2*tau1) + 3*(RC.b31*pow(Ax_li,3.0) - RC.b32*Ax_li*pow(Az_li,2.0))*cos(3*tau1)   );
    yv_li[5] = dtau1dt*( -dm*Az_li*sin(tau1) -2*dm*RC.d21*Ax_li*Az_li*sin(2*tau1) - 3*dm*(RC.d32*Az_li*pow(Ax_li,2.0) - RC.d31*pow(Az_li,3.0))*sin(3*tau1)  );


    //Back into COM frame
    switch(halo.Li.number) //special case of the x dimension
    {
    case 1:
        yv_li[0] = yv_li[0]*gamma_i - mu + 1 - gamma_i;
        break;
    case 2:
        yv_li[0] = yv_li[0]*gamma_i - mu + 1 + gamma_i;
        break;
    case 3:
        yv_li[0] = yv_li[0]*gamma_i - mu - gamma_i;
        break;
    default: //L2
        yv_li[0] = yv_li[0]*gamma_i - mu + 1 + gamma_i;
        break;
    }
    for(i=1; i<6; i++)yv_li[i]*=gamma_i; //other dimensions
}


void apply_halo_orbit_third_order(Orbit *orbit)
{
    //Init
    gsl_matrix *STM0 = gsl_matrix_alloc(6,6);
    gsl_matrix_set_identity (STM0);

    //Third order approximation of the state
    third_order_halo_orbit(*orbit, 0.0, orbit->y0);

    //STM(0) is concatenated to the state @t=0
    custom_matrixToVector(orbit->y0, STM0, 6, 6, 6);
    gsl_matrix_free(STM0);

    //Updating orbital period
    orbit->T = third_order_period(*orbit);

    //Half-period value
    third_order_halo_orbit(*orbit, 0.5*orbit->T, orbit->yhalf);

    //Status
    orbit->status = HALO_APPROXIMATED;
}


double third_order_period(Orbit orbit)
{
double gamma_i = orbit.Li.gamma_i;


//Richardson coefficients
    RichardsonCoeffs RC;
    richardson_coefficients(orbit.cr3bp, orbit.Li.number, &RC);

//X and Z-amplitude in Li-frame
    double Az_li = orbit.Az_R/gamma_i;
    double Ax_li = sqrt(-(RC.l2*pow(Az_li,2.0) + RC.Delta)/RC.l1);

//nu parameter
    double nu = 1+RC.s1*pow(Ax_li,2.0)+RC.s2*pow(Az_li,2.0);

//Third order estimate of the orbital period
    return (2*M_PI)/(RC.omega_p*nu);

}


void richardson_coefficients(CR3BP cr3bp, int pointNumber, RichardsonCoeffs *RC)
{
    //---------------------------------------------------------------------------------------------------------------------------
    //Init
    //---------------------------------------------------------------------------------------------------------------------------
    double gamma_i;
    double mu = cr3bp.mu;

    switch(pointNumber)
    {
    case 1:
        gamma_i = cr3bp.l1.gamma_i;
        break;
    case 2:
        gamma_i = cr3bp.l2.gamma_i;
        break;
    case 3:
        gamma_i = cr3bp.l3.gamma_i;
        break;
    default:
        gamma_i = cr3bp.l1.gamma_i;
        printf("gamma_i is initialized to its default value (l1-m2)");
        break;
    }


    //Cn coefficients for the corresponding li point
    double c2 = cn(mu, gamma_i, pointNumber, 2);
    double c3 = cn(mu, gamma_i, pointNumber, 3);
    double c4 = cn(mu, gamma_i, pointNumber, 4);

    RC->c2 = c2;
    RC->c3 = c3;
    RC->c4 = c4;

    //---------------------------------------------------------------------------------------------------------------------------
    //Computation of characteristic constants of the periodic motion
    //---------------------------------------------------------------------------------------------------------------------------
    //lambda
    double lambda;
    lambda = sqrt(0.5*(2 - c2 + sqrt(9 * pow(c2,2.0) - 8*c2)));
    RC->lambda = lambda;
    //omega_p
    RC->omega_p = sqrt(0.5*(2 - c2 + sqrt(9 * pow(c2,2.0) - 8*c2)));  //omega_p is set equal to lambda, see Koon et al. & Gomez et al. for details
    //omega_v
    RC->omega_v = sqrt(c2);
    //kappa
    double kappa = (pow(lambda,2.0)+1+2*c2)/(2.0*lambda);
    RC->kappa = kappa;
    //Delta
    RC->Delta = pow(lambda,2.0)-c2;

    //Coeffs
    double d1 = (3*pow(lambda,2.0))/kappa*( kappa*(6*pow(lambda,2.0) - 1) - 2*lambda);
    double d2 = (8*pow(lambda,2.0))/kappa*( kappa*(11*pow(lambda,2.0) - 1) - 2*lambda);

    RC->d1 = d1;
    RC->d2 = d2;

    RC->a21 = (3*c3*(pow(kappa,2.0) - 2))/(4*(1 + 2*c2));
    RC->a22 = (3*c3)/(4*(1 + 2*c2));
    RC->a23 = - (3*c3*lambda)/(4*kappa*d1)*(3*pow(kappa,3.0)*lambda - 6*kappa*(kappa-lambda) + 4);
    RC->a24 = - (3*c3*lambda)/(4*kappa*d1)*(2 + 3*kappa*lambda);

    RC->b21 = - (3*c3*lambda)/(2*d1)*(3*kappa*lambda - 4);
    RC->b22 = (3*c3*lambda)/(d1);

    RC->d21 = - c3/(2*pow(lambda,2.0));

    RC->a31 =  - (9*lambda)/(4.0*d2)*(4*c3*(kappa*(RC->a23) - (RC->b21)) + kappa*c4*(4 + pow(kappa,2.0))) + (9*pow(lambda,2.0) + 1 - c2)/(2.0*d2)*(3*c3*(2*RC->a23 - kappa*RC->b21) + c4*(2+3*pow(kappa,2.0)));
    RC->a32 = - (9*lambda)/(4.0*d2)*(4*c3*(kappa*RC->a24 - RC->b22) + kappa*c4) - 3/(2*d2)*(9*pow(lambda,2.0) + 1 - c2)*(c3*(kappa*RC->b22 + RC->d21 - 2*RC->a24) - c4);

    RC->b31 = 3/(8*d2) * ( 8*lambda*(3*c3*(kappa*RC->b21 - 2*RC->a23) - c4*(2 + 3*pow(kappa,2.0))) + (9*pow(lambda,2.0) + 1 + 2*c2)*(4*c3*(kappa*RC->a23 - RC->b21) + kappa*c4*(4 + pow(kappa,2.0))) );
    RC->b32 =  1.0/d2*(9*lambda*(c3*(kappa*RC->b22 + RC->d21 - 2*RC->a24) - c4) + 3/8.0*(9*pow(lambda,2.0) + 1 + 2*c2)*(4*c3*(kappa*RC->a24 - RC->b22) + kappa*c4));

    RC->d31 = 3/(64*pow(lambda,2.0))*(4*c3*RC->a24 + c4);
    RC->d32 = 3/(64*pow(lambda,2.0))*(4*c3*(RC->a23 - RC->d21) + c4*(4 + pow(kappa,2.0)));

    RC->s1 = pow(2*lambda*(lambda*(1 + pow(kappa,2.0)) - 2*kappa),-1.0) * (3/2.0*c3*(2*RC->a21*(pow(kappa,2.0) - 2) -  RC->a23*(pow(kappa,2.0) + 2) - 2*kappa*RC->b21)- 3/8.0*c4*(3*pow(kappa,4.0) - 8*pow(kappa,2.0) + 8));
    RC->s2 = pow(2*lambda*(lambda*(1 + pow(kappa,2.0)) - 2*kappa),-1.0) * (3/2.0*c3*(2*RC->a22*(pow(kappa,2.0) - 2) + RC->a24*(pow(kappa,2.0) + 2) + 2*kappa*RC->b22 + 5*RC->d21) + 3/8.0*c4*(12 - pow(kappa,2.0)));


    RC->l1 = - 3/2.0 * c3*(2*RC->a21 + RC->a23 +5*RC->d21) - 3/8.0*c4*(12 - pow(kappa,2.0)) + 2*pow(lambda,2.0)*RC->s1;
    RC->l2 =   3/2.0 * c3*(RC->a24 - 2*RC->a22) + 9/8.0*c4 + 2*pow(lambda,2.0)*RC->s2;

}


void print_richardson_coefficients(CR3BP cr3bp, int pointNumber)
{
    //Richardson coefficients
    RichardsonCoeffs RC;
    richardson_coefficients(cr3bp, pointNumber, &RC);

    printf("------------------------------------------\n");
    printf("Richardson Coefficients for L%i\n", pointNumber);
    printf("------------------------------------------\n");


    printf("lamb: %+.5e\n", RC.lambda);
    printf("kap:  %+.5e\n", RC.kappa);

    printf("D:    %+.5e\n", RC.Delta);

    printf("c2:   %+.5e\n", RC.c2);
    printf("c3:   %+.5e\n", RC.c3);
    printf("c4:   %+.5e\n", RC.c4);

    printf("s1:   %+.5e\n", RC.s1);
    printf("s2:   %+.5e\n", RC.s2);

    printf("l1:   %+.5e\n", RC.l1);
    printf("l2:   %+.5e\n", RC.l2);

    printf("d1:   %+.5e\n", RC.d1);
    printf("d2:   %+.5e\n", RC.d2);

    printf("a21:  %+.5e\n", RC.a21);
    printf("a22:  %+.5e\n", RC.a22);
    printf("a23:  %+.5e\n", RC.a23);
    printf("a24:  %+.5e\n", RC.a24);

    printf("a31:  %+.5e\n", RC.a31);
    printf("a32:  %+.5e\n", RC.a32);

    printf("b21:  %+.5e\n", RC.b21);
    printf("b22:  %+.5e\n", RC.b22);

    printf("b31:  %+.5e\n", RC.b31);
    printf("b32:  %+.5e\n", RC.b32);

    printf("d21:  %+.5e\n", RC.d21);
    printf("d31:  %+.5e\n", RC.d31);
    printf("d32:  %+.5e\n", RC.d32);

    printf("------------------------------------------\n");
}


double cn(double mu, double gamma_i, int pointNumber, int n)
{
    switch(pointNumber)
    {
    case 1:
        return pow(gamma_i,-3.0) * ( mu + pow(-1,n)*(1-mu)*pow(gamma_i, n+1)/pow(1 - gamma_i, n+1) );
    case 2:
        return pow(gamma_i,-3.0) * ( pow(-1.0,n)*mu + pow(-1,n)*(1-mu)*pow(gamma_i, n+1)/pow(1 + gamma_i, n+1) );
        break;
    case 3:
        return pow(gamma_i,-3.0) * ( 1 - mu + mu*pow(gamma_i,n+1)/pow(1 + gamma_i, n+1) );
        break;
    }
    return 0; //never here
}
