#include <math.h>
#include <gsl/gsl_matrix.h>
#include "custom.h"
/**
* \brief custom routine to transform vector to matrix with a given shift in the initial vector
**/
void custom_vectorToMatrix(gsl_matrix *m, double y[], int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            gsl_matrix_set(m, i-1, j-1, y[shift+k-1]);
        }
}

/**
* \brief custom routine to transform matrix to vector with a given shift in the final vector
**/
void custom_matrixToVector(double y[], gsl_matrix *m, int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            y[shift+k-1] = gsl_matrix_get(m, i-1, j-1);
        }
}
