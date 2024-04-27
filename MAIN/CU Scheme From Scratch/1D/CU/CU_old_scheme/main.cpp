#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>
#include<sstream>
#include<fstream>

#include "Constants.h"
#include "Utility.h"
#include "OutputResult.h"
#include "InitializeConservedVariables.h"
#include "ExtendCells.h"
#include "CUNumericalFluxOld.h"
#include "SSPRK.h"

using namespace std;

namespace cts = Constants;
namespace utl = Utility;
namespace opr = OutputResult;
namespace exc = ExtendCells;
namespace cuf = CUNumericalFLux;
namespace icv = InitializeConservedVariables;
namespace srk = SSPRK;

int main(){
    /*
        ! Export run environment for later use
    */
    utl::export_env();

    /*
        ! Definitions for the required variables
    */
    double T = cts::time.first;
    double dt;

    //! Define computational grid and the conserved variable vectors
    // computational grid
    vd grid(cts::n+2, 0);

    // vector to store all 3 conserved variables [density, momentum, energy]
    vvd U(3, vd(cts::n+2, 0));

    /*
        ! Simulation steps start from now
    */
    //! initialize computational grid
    for(int i=1 ; i<=cts::n ; i++){
        grid[i] = cts::dom.first + cts::dx*( i - 0.5 );
    }

    // output the computational grid
    opr::write_grids(cts::MODE, grid);

    /*
        ! Initialize the conserved variables and extend the cells 
        ! based on appropriate boundary conditions
    */
    icv::initialize_conserved_variables(U, grid);

    // output initial conditions
    opr::write_conserved_variables(cts::MODE, U[0]);

    /*
        ! Now we can start the iterations
    */    

    while( T < cts::time.second ){
        cout << "T=" << T << " | dt=" << dt << endl;
        srk::update_conserved_variables(U, dt, T);
        T += dt;
    }

    cout << "Final Time (T) = " << T << endl;

    // output final conditions
    opr::write_conserved_variables(cts::MODE, U[0]);

    return 0;
}
