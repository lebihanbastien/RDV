#ifndef CUSTOM_H_INCLUDED
#define CUSTOM_H_INCLUDED

#include <gsl/gsl_matrix.h>


void custom_vectorToMatrix(gsl_matrix *m, double y[], int rows, int columns, int shift);
void custom_matrixToVector(double y[], gsl_matrix *m, int rows, int columns, int shift);

#endif // CUSTOM_H_INCLUDED
