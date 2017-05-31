#include "vf.h"



//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations)
int qbfbp_vfn_novar(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @ t
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, qbp->nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, qbp->nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, qbp->nf, qbp->cs.pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    return 0;
}

//Update the normalized vector field of the state
int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma)
{
    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = alpha[0]*y[3] + alpha[1]*y[0] + alpha[2]*y[1];
    //f[1] = y'
    //------------------------------------
    f[1] = alpha[0]*y[4] + alpha[1]*y[1]-alpha[2]*y[0];
    //f[2] = z'
    //------------------------------------
    f[2] = alpha[0]*y[5] + alpha[1]*y[2];
    //f[3] = px'
    //------------------------------------
    f[3] = - alpha[1]*y[3] + alpha[2]*y[4] + alpha[14]*y[0]
           + alpha[12];
    if(me != 0) f[3] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[0] - pe[0]);
    if(mm != 0) f[3] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[0] - pm[0]);
    if(ms != 0) f[3] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[0] - ps[0]);

    //f[4] = py'
    //------------------------------------
    f[4] = - alpha[1]*y[4] - alpha[2]*y[3] + alpha[14]*y[1]
           + alpha[13];
    if(me != 0) f[4] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[1] - pe[1]);
    if(mm != 0) f[4] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[1] - pm[1]);
    if(ms != 0) f[4] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[1] - ps[1]);

    //f[5] = pz'
    //------------------------------------
    f[5] = - alpha[1]*y[5] + alpha[14]*y[2];
    if(me != 0) f[5] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[2] - pe[2]);
    if(mm != 0) f[5] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[2] - pm[2]);
    if(ms != 0) f[5] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[2] - ps[2]);

    return 0;
}
