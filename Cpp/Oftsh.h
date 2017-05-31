#ifndef OFTSH_H
#define OFTSH_H


class Oftsh
{
    private:

    int nv;            ///number of variables
    int order;         ///order

    int cnv;           ///number of variables of the coefficients
    int corder;        ///order of the coefficients

    T *coefs;          ///coefficients in the form of Fourier series
};

#endif // OFTSH_H
