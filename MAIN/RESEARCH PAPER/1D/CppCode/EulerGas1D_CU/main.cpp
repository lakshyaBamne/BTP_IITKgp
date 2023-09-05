/*
    1D Central Upwind Scheme to solve Riemann Problems
*/

#include<iostream>
#include<vector>
#include<utility>
#include<string>

#include "InitializeRiemannProblem.h"
#include "CUNumericalFlux.h"
#include "SSPRK.h"
#include "OutputResult.h"
#include "EulerForwardDifference.h"
#include "PrimitiveVariables.h"

#define ll long long int

using namespace std;

namespace IRP = InitRiemannProblem;
namespace CUF = CUNumericalFlux;
namespace SRK = SSPRK;
namespace OPF = OutputResult;
namespace EFD = EulerForwardDifference;
namespace PRV = PrimitiveVariables;

int main(){
    
    // Variables required to take initial conditions as input from the user
    pair<double,double> domain; // { x_start, x_end }
    pair<double,double> time; // { t_start, t_end }
    ll Nx; // number of grid points in the computational domain
    pair<string,string> initial_conditions; // riemann problem and boundary condition
    
    cout << "-----------------INPUT INITIAL CONDITIONS--------------" << endl;
    
    cout << "Domain start : ";
    cin >> domain.first;

    cout << "Domain end : ";
    cin >> domain.second;

    cout << "Grid points in domain : ";
    cin >> Nx;

    cout << "Initial time : ";
    cin >> time.first;

    cout << "Final time : ";
    cin >> time.second;

    cout << "Enter the Riemann Problem : ";
    cin >> initial_conditions.first;

    cout << "Enter the Boundary Conditions : ";
    cin >> initial_conditions.second;

    double dt;
    double t=0; // initial time
    double dx = (domain.second - domain.first)/Nx;

    // STEP 1 Initialize the Conserved variable vectors using the Initial conditions
    // -> Finite Volume Grid is created
    // -> Conserved variables are initialized according to the Riemann Problem
    // -> Ghost values are added according to the given boundary conditions
    vector<vector<double>> cons_vars = IRP::get_conserved_variables(initial_conditions, domain, Nx);

    // //! Writing output to a file for plotting later
    // OPF::write_vector(cons_vars[0], "Density2.txt");
    // OPF::write_vector(cons_vars[1], "Momentum2.txt");
    // OPF::write_vector(cons_vars[2], "Energy2.txt");

    // // convert the conserved variables to primitive variables for plotting
    // vector< vector<double> > prim_vars = PRV::get_primitive_variables(cons_vars);

    // OPF::write_vector(prim_vars[0], "Velocity2.txt");
    // OPF::write_vector(prim_vars[1], "Pressure2.txt");

    // Now we can start the iterations
    while( t < time.second ){
        // Log
        cout << "t = " << t << " | dt = " << dt << endl;

        //! Writing output to a file for plotting later
        OPF::write_vector(cons_vars[0], "Density2.txt");
        OPF::write_vector(cons_vars[1], "Momentum2.txt");
        OPF::write_vector(cons_vars[2], "Energy2.txt");

        // convert the conserved variables to primitive variables for plotting
        vector< vector<double> > prim_vars = PRV::get_primitive_variables(cons_vars);

        OPF::write_vector(prim_vars[0], "Velocity2.txt");
        OPF::write_vector(prim_vars[1], "Pressure2.txt");

        vector<vector<double>> cu_flux = CUF::get_cu_flux(cons_vars, initial_conditions, dt, dx, t, time);

        double LAMBDA = dt / dx;
        
        t = t+dt;

        vector<vector<double>> cons_vars_next = SRK::get_next_cons_vars(cons_vars, cu_flux, LAMBDA, initial_conditions, dt, dx, t, time);    
        // vector<vector<double>> cons_vars_next = EFD::get_next_cons_vars(cons_vars, cu_flux, LAMBDA, initial_conditions, dt, dx, t, time);

        // copy the new values in the old vector to be used in the next iteration
        for(int i=0 ; i<cons_vars[0].size() ; i++){
            cons_vars[0][i] = cons_vars_next[0][i];
            cons_vars[1][i] = cons_vars_next[1][i];
            cons_vars[2][i] = cons_vars_next[2][i];
        }
    }

    // //! Writing output to a file for plotting later (only the final values)
    // OPF::write_vector(cons_vars[0], "Density2.txt");
    // OPF::write_vector(cons_vars[1], "Momentum2.txt");
    // OPF::write_vector(cons_vars[2], "Energy2.txt");

    // prim_vars = PRV::get_primitive_variables(cons_vars);

    // OPF::write_vector(prim_vars[0], "Velocity2.txt");
    // OPF::write_vector(prim_vars[1], "Pressure2.txt");

    return 0;
}