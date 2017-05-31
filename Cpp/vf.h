#ifndef VF_H_INCLUDED
#define VF_H_INCLUDED

#include "env.h"
#include "pmcoc.h"

//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations)
int qbfbp_vfn_novar(double t, const double y[], double f[], void *params_void);

int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma);

#endif // VF_H_INCLUDED
