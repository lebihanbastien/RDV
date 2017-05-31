#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

/**
 * \file nrutils.h
 * \brief Useful routines from NR in C.
 */


#include <complex.h>
/*
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
*/

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror(char error_text[]);
float *fvector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
double complex *dcvector(long nl, long nh);
float **cmatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
double complex **dcmatrix(long nrl, long nrh, long ncl, long nch);
int **icmatrix(long nrl, long nrh, long ncl, long nch);
float **subcmatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_cmatrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_fvector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_dcvector(double complex *v, long nl, long nh);
void free_cmatrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dcmatrix(double complex **m, long nrl, long nrh, long ncl, long nch);
void free_icmatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_subcmatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_cmatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

#else /* ANSI */
/* traditional - K&R */

void nrerror();
float *fvector();
float **cmatrix();
float **subcmatrix();
float **convert_cmatrix();
float ***f3tensor();
double ***f3tensor();
double *dvector();
double complex *dcvector();
double **dmatrix();
double complex **dcmatrix();
int *ivector();
int **icmatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_fvector();
void free_dvector();
void free_dcvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_cmatrix();
void free_subcmatrix();
void free_convert_cmatrix();
void free_dmatrix();
void free_dcmatrix();
void free_icmatrix();
void free_f3tensor();
void free_d3tensor();

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */
