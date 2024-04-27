#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<cmath>

#include "Constants.h"
#include "ExtendCells.h"

using namespace std;

namespace cts = Constants;
namespace exc = ExtendCells;

namespace InitializeConservedVariables{
    void initialize_conserved_variables(vvd& U, vd grid);
}

void InitializeConservedVariables::initialize_conserved_variables(vvd& U, vd grid){

    double u, p;

    if( cts::PROBLEM == "MCW" ){

        for(int i=1 ; i<=cts::n ; i++){
            if( grid[i] <= 0.3 ){
                U[0][i] = 1.4;
                u = 0.,1;
                p = 1;
            }
            else{
                U[0][i] = 1.0;
                u = 0.1;
                p = 1;
            }

            U[1][i] = U[0][i]*u; 
            U[2][i] = p/(cts::GAMMA-1) + 0.5*U[0][i]*pow(u,2);
        }

    }
    else if( cts::PROBLEM == "SCW" ){



    }
    else if( cts::PROBLEM == "BLW" ){

        for(int i=1 ; i<=cts::n ; i++){

            if( grid[i] < 0.1 ){
                U[0][i] = 1.0;
                u = 0;
                p = 1000;
            }
            else if( grid[i]>=0.1 && grid[i]<=0.9 ){
                U[0][i] = 1.0;
                u = 0;
                p = 0.01;
            }
            else{
                U[0][i] = 1.0;
                u = 0;
                p = 100;
            }

            U[1][i] = U[0][i]*u;
            U[2][i] = p/(cts::GAMMA-1) + 0.5*U[0][i]*pow(u,2);
        }

    }
    else if( cts::PROBLEM == "SOS" ){



    }
    else if( cts::PROBLEM == "LAX" ){



    }
    else{
        cout << "---ERROR--- Please enter correct problem ---" << endl;
    }

    //! extend cells using the appropriate boundary conditions for the problem
    exc::extend_conserved_variables(U);

}
