#include<iostream>
#include<vector>
#include<utility>
#include<string>

#include "Parameters.h"
#include "Constants.h"
#include "Utility.h"
#include "Initialize.h"
#include "Output.h"
#include "ExtendCells.h"
#include "CU1D.h"
#include "SSPRK.h"

using namespace std;

namespace UTL = Utility;
namespace CTS = Constants;
namespace PRM = Parameters;
namespace INI = Initialize;
namespace OUT = Output;
namespace EXC = ExtendCells;
namespace CUS = CentralUpwindScheme;
namespace SRK = SSPRK;

#define ll long long int
#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define pdd pair<double,double>
#define pss pair<string,string>
#define pii pair<ll,ll>

int main(){

    /*
        Variable declarations used in the simualtion
    */
    pdd domain;
    pdd time;
    ll Nx;

    double dt;

    // Initialize the Domain variables
    INI::InitializeVariables(domain, time, Nx);

    double t = time.first;

    // dx value
    double dx = (domain.second-domain.first)/Nx;

    cout << "-----------Initialized Domain Variables----------" << endl;
    cout << "Domain : [ " << domain.first << " , " << domain.second << " ]" << endl;
    cout << "Time : [ " << time.first << " , " << time.second << " ]" << endl;
    cout << "Number of Finite Volume Cells : " << Nx << endl;

    // Initialize the Computational Grid
    vd grid = INI::MakeGrid(domain, Nx);

    cout << "----------Initialized Computational Grid-----------" << endl;
    cout << "x = { ";
    for(auto i : grid){
        cout << i << ", ";
    }
    cout << " }" << endl;

    // Write initial conditions to output file for later use
    OUT::write_computational_grid(grid);

    // Initialize the Riemann Problem Variables
    vvd cons_vars = INI::InitializeRiemannProblem(grid);

    cout << "---------Initialized Riemann Conserved Variables----------" << endl;

    // Now we are ready to start the Central Upwind Scheme

    cout << "--- STARTING SIMULATIONS ---" << endl;

    switch (PRM::SIMTYPE){
        case 1: // PLOT

            OUT::write_conserved_variables(cons_vars);

            while( t < time.second ){
                cout << "t = " << t << " | dt = " << dt << endl;
                SRK::update_conserved_variables(cons_vars, dt, dx);
                t = t + dt;
            }        

            OUT::write_conserved_variables(cons_vars);

            break;
        
        case 2: // ANIMATE

            while( t < time.second ){
                cout << "t = " << t << " | dt = " << dt << endl;
                
                OUT::write_conserved_variables(cons_vars);
                SRK::update_conserved_variables(cons_vars, dt, dx);
                t = t + dt;
            }

            OUT::write_conserved_variables(cons_vars);

            break;

        default:
            cout << "---ERROR--- Please enter correct Simulation type ---" << endl;
            break;
    }

    return 0;
}

