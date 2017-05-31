#ifndef MANIFOLD_H_INCLUDED
#define MANIFOLD_H_INCLUDED


typedef struct Manifold_branch Manifold_branch;
struct Manifold_branch
{
    Orbit orbit;                    //Original orbit
    int stability;                  //Stability type
    int way;                        //Interior or exterior
    struct value_params val_par;    //Event parameters
    double final_time;              //final time after computation
    double theta;                   //Position of departure on the orbit
    gsl_vector *direction;          //Direction of departure
    double initial_position[6];     //Initial position
    double final_position[6];       //Final position
    double** events_mat;            //Pointer towards the stored position of events;
    struct value_function fvalue;   //fvalue for event detection
};

void init_manifold_branch(Manifold_branch *branch, Orbit orbit, int stability, int way, struct value_params val_par);
void free_manifold_branche(Manifold_branch branch);
int manifold_branch_computation(Manifold_branch *branch,
                                double theta,
                                double tf,
                                custom_ode_structure *ode_s,
                                custom_ode_structure *ode_s_6);
void manifold_branch_plot(Manifold_branch *branch, gnuplot_ctrl *h1, int type, int points, int scale, custom_ode_structure *ode_s_6);
#endif // MANIFOLD_H_INCLUDED
