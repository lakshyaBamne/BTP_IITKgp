#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<unordered_map>

#include "Constants.h"
#include "OutputResult.h"
#include "TimeIntegration.h"
#include "InitializeComputationalGrid.h"
#include "InitializeConservedVariables.h"

using namespace std;

namespace cts = Constants;
namespace opr = OutputResult;
namespace srk = SSPRK;
namespace efd = EFD;
namespace icg = InitializeComputationalGrid;
namespace icv = InitializeConservedVariables;

int main(){
    cout << "---RUNNING 1D SCHEMES---" << endl;
    
    string mode;
    
    cout << "-> [PLOT-ONE] : only plot a single instance of the problem" << endl;
    cout << "-> [PLOT-REF] : plot two instances of the problem (one being reference)" << endl;
    cout << "Enter mode: ";
    cin >> mode;

    if( mode == "PLOT-ONE" ){
        opr::export_env("ENV/env.txt");
    }
    else if( mode == "PLOT-REF" ){
        opr::export_env("ENV/env1.txt", "ENV/env2.txt");
    }
    else{
        cout << "---ERROR--- Please enter correct mode to run the scheme ---" << endl;
    }
    
    // start the simulations
    double dt;
    double T = cts::time.first;

    // Initialize computational grid and conserved variables
    vd grid = icg::make_grid();

    opr::output_vector(grid, "GRIDS/"+cts::PROBLEM+".txt");

    vvd U = icv::initialize_conserved_variables(grid);

    // output initial conditions
    string path = "RESULTS/" + cts::PROBLEM + "-" + cts::FLUX + ".txt";
    opr::output_vector(U[0], path);

    while( T < cts::time.second ){
        cout << "T=" << T << " | dt=" << dt << endl;

        if( cts::TI == "SSPRK" ){
            srk::update_conserved_variables(U, dt, T);
        }
        else if( cts::TI == "EFD" ){
            efd::update_conserved_variables(U, dt, T);
        }
        else{
            cout << "---ERROR--- Please enter correct Time integration scheme ---" << endl;
        }
        
        T += dt;
    }

    // output final conditions
    opr::output_vector(U[0], path);

    return 0;
}