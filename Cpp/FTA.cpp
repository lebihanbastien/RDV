#include "FTA.h"

#include <stdio.h>
#include <stdlib.h>

/**
 * \file ftda.cpp
 * \brief Basic operations for algebraic manipulation (exponents in reverse lexicographic order, number of coefficient in homogeneous polynomials...).
 * \author BLB, Angel Jorba
 * \date Feb 2016
 * \version 1.0
 *
 *      Based on Jorba 1999.
 */

using namespace std;

int  FTA::nor=0;       //static maximum degree
int  FTA::nov=0;       //static number of variables
long int **FTA::psi;   //cointains psi(i,j) for i=2..nov.


//-----------------------------------------------------------------------
//Class routines
//-----------------------------------------------------------------------
/**
 *   \brief Initializes the table psi(i,j), which contains the number of monomials of degree j with i variables.
 *
 *   parameters:
 *   nr: maximum degree we are going to work with. it can not be greater than 63.
 *   nv: number of variables
 *
 *   returned value: number of kbytes allocated by the internal tables.
 *   Based on a routine by Angel Jorba, 1999.
 **/
int  FTA::init(int nv, int nr)
{
    int i,j,l;
    unsigned long int mem; /* mem: to count the amount of memory used */

    if (nor != 0) /* this means that FTA is already initialized */
    {
        if (nr <= nor)
        {
            puts("**************************************************");
            printf("%s warning message:\n"," init");
            printf("%s is already initialized to degree %d\n"," init",nor);
            printf("now you want to initialize it again to the same degree.\n");
            printf("action taken: this call to %s is simply ignored.\n"," init");
            puts("**************************************************");
            return(0);
        }
        printf("%s error message:\n"," init");
        printf("%s is already initialized to degree %d\n"," init",nor);
        printf("and now you want to initialize it to degree %d.\n",nr);
        printf("you must call routine %s first.\n"," free");
        printf("action taken: program aborted\n");
        exit(1); /* this will stop the program */
    }
    nor=nr;
    nov=nv;
    mem=(nv-1)*sizeof(int*); /* this is to count the number of kb. allocated */
    psi= (long  int**)malloc((nv)*sizeof(long int*));
    if (psi == NULL)
    {
        puts("FTA error. no memory (1).");
        exit(1);
    }
    psi -= sizeof(long  int*);
    for (i=1; i<=nv; i++)
    {
        psi[i]=(long int*)malloc((nr+1)*sizeof(long  int));
        if (psi[i] == NULL)
        {
            puts("FTA error. no memory (2).");
            exit(1);
        }
        mem += (nr+1)*sizeof(int);
    }
    for (j=0; j<=nr; j++) psi[1][j]=1;

    if(nv > 1)
    {
        for (j=0; j<=nr; j++) psi[2][j]=j+1;
        for (i=3; i<=nv; i++)
        {
            for (j=0; j<=nr; j++)
            {
                psi[i][j]=0;
                for (l=0; l<=j; l++) psi[i][j] += psi[i-1][l];
            }
        }
    }
    mem += (nr+1)*sizeof(int*);
    mem /= 1024;

    return(mem);
}

/**
 *  \brief Frees the space allocated by FTA::init.
 **/
void  FTA::free()
{
    int i;
    if (nor == 0)
    {
        puts("**************************************************");
        printf("%s warning message:\n"," free");
        printf("no memory to free\n");
        printf("action taken: this call is simply ignored.\n");
        puts("**************************************************");

        return;
    }
    cout << "nov: " << nov << endl;
    for (i=1; i<=nov; i++)
    {
        cout << "i: " << i << endl;
        delete psi[i];
    }
    psi += sizeof(long  int*);
    delete psi;
    nor=0;
    nov=0;
    return;
}

/**
 *   \brief Returns the number of monomials of degree nr with nv variables, making use of the table FTA::psi.
 *          For now, uses the routine binomial rather than the matrix psi. Hence, FTA::int does NOT need to be called.
 *
 *  \param  nv: number of variables
 *  \param  nr: order we are interested in (input).
 *  \return value: number of monomials of order no.
 *
 **/
long int FTA::nmon(int nv, int nr)
{
//    if (nr > nor)
//    {
//        printf("nmon: error, the requested degree %i is greater than nor=%i.\n", nr, nor );
//        exit(1);
//    }
//    if (nv > nov)
//    {
//        puts("nmon: error, the requested numver of variable is greater than nov.");
//        exit(1);
//    }

    if(nv <= 1)
        return 1;
    else
    return(binomial(nv+nr-1, nv-1));
    //return(psi[nv][nr]);

}

/**
 *  \brief given a multiindex k, this routine computes the next one
 *  according to the (reverse) lexicographic order.
 *
 *  \param k: array of nv components containing the multiindex. it is overwritten on exit (input and output).
 **/
void  FTA::prxkt(int k[], int nv)
{
    if(nv == 0)
    {
        k[0]++;
        return;
    }
    else
    {
        int i;
        if (k[0] != 0)
        {
            k[0]--;
            k[1]++;
            return;
        }
        for(i=1; i<nv-1; i++)
        {
            if (k[i] != 0)
            {
                k[0]=k[i]-1;
                k[i]=0;
                k[i+1]++;
                return;
            }
        }
    }
    puts("prxkt error 1.");
    exit(1);
}

/**
 *   \brief Returns the number of product operation in a taylor serie product
 **/
long int pdk(int nv, int n)
{
    long int result = 0;
    int k,i;
    for(k=0; k<=n; k++)
    {
        for(i=0; i<=k; i++) result+= FTA::nmon(nv,i)* FTA::nmon(nv,k-i); //result+=binomial(nv+i-1, nv-1)*binomial(nv+k-i-1, nv-1); //
    }
    return result;
}


//-----------------------------------------------------------------------
//Binomial coefficients
//-----------------------------------------------------------------------
/**
 *  \brief Inner routine for the computation of the binomial (n k) in the routine binomial.
 **/
unsigned long gcd_ui(unsigned long x, unsigned long y)
{
    unsigned long t;
    if (y < x)
    {
        t = x;
        x = y;
        y = t;
    }
    while (y > 0)
    {
        t = y;
        y = x % y;
        x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
    }
    return x;
}

/**
 *  \brief Returns the binomial (n k)
 **/
unsigned long binomial(unsigned long n, unsigned long k)
{
    unsigned long d, g, r = 1;
    if (k == 0) return 1;
    if (k == 1) return n;
    if (k >= n) return (k == n);
    if (k > n/2) k = n-k;
    for (d = 1; d <= k; d++)
    {
        if (r >= ULONG_MAX/n)    /* Possible overflow */
        {
            unsigned long nr, dr;  /* reduced numerator / denominator */
            g = gcd_ui(n, d);
            nr = n/g;
            dr = d/g;
            g = gcd_ui(r, dr);
            r = r/g;
            dr = dr/g;
            if (r >= ULONG_MAX/nr) return 0;  /* Unavoidable overflow */
            r *= nr;
            r /= dr;
            n--;
        }
        else
        {
            r *= n--;
            r /= d;
        }
    }
    return r;
}


