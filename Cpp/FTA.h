#ifndef FTA_H
#define FTA_H

#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

class FTA
{
    ///static maximum degree
    static int nor;
    ///static number of variables
    static int nov;
    ///contains psi(i,j) for i=2..nov, which is the number of monomials of degree j with i variables.
    static long  int **psi;

public:

    /**
     *   \brief Initializes the table psi(i,j), which contains the number of homogeneous polynomials of order j with i variables.
     *
     *   parameters:
     *   nr: maximum degree we are going to work with. it can not be greater than 63.
     *   nv: number of variables
     *
     *   returned value: number of kbytes allocated by the internal tables.
     *   Based on a routine by Angel Jorba, 1999.
     **/
    static int init(int nv, int nr);
    /**
     *  \brief Frees the space allocated by FTDA::init.
     **/
    static void free();
    /**
     *   \brief Returns the number of monomials of degree nr with nv variables, making use of the table FTDA::psi.
     *
     *  parameters:
     *  nv: number of variables
     *  nr: order we are interested in (input).
     *  returned value: number of monomials of order no.
     **/
    static long int nmon(int nv, int nr);
    /**
     *  \brief given a multiindex k, this routine computes the next one
     *  according to the (reverse) lexicographic order.
     *
     *  parameters:
     *  k: array of nv components containing the multiindex. it is overwritten on exit (input and output).
     **/
    static void prxkt(int k[], int nv);
};

//-----------------------------------------------------------------------
//Binomial coefficients
//-----------------------------------------------------------------------

unsigned long binomial(unsigned long n, unsigned long k);
unsigned long gcd_ui(unsigned long x, unsigned long y);

#endif // FTA_H
