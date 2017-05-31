#ifndef USER_PARAMETERS_H_INCLUDED
#define USER_PARAMETERS_H_INCLUDED

//Differential correction structure
typedef struct diff_corr_parameters diff_corr_parameters;
struct diff_corr_parameters
{
    double eps_diff;     //for differential correction
    int diff_corr_type;  //type of differential correction
};

//Ode structure
typedef struct ode_parameters ode_parameters;
struct ode_parameters
{
    double eps_abs;     //Abs precision
    double eps_rel;     //Rel precision
};

//Global user structure
typedef struct user_parameters user_parameters;
struct user_parameters
{
    //Differential correction structure
    diff_corr_parameters diff_corr;
    ode_parameters ode;
};



#endif // USER_PARAMETERS_H_INCLUDED
