/*
    ! Namespace to initialize the conserved variables based on the problem
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace InitializeConservedVariables{
    vvd initialize_conserved_variables(const vd& grid);
}

vvd InitializeConservedVariables::initialize_conserved_variables(const vd& grid){
    vvd U(3, vd(cts::n+2, 0));
    
    double u, p;

    if( cts::PROBLEM == "MCW" ){

        for(int i=1 ; i<=cts::n ; i++){
            if( grid[i] <= 0.3 ){
                U[0][i] = 1.4;
                u = 0.1;
                p = 1.0;
            }
            else{
                U[0][i] = 1.0;
                u = 0.1;
                p = 1.0;
            }

            U[1][i] = U[0][i]*u;
            U[2][i] = ( p/(cts::GAMMA-1) ) + ( U[0][i]*pow(u,2)/2 );
        }

    }
    else if( cts::PROBLEM == "SCW" ){

        for(int i=1 ; i<=cts::n ; i++){
            if( grid[i] <= 0.8 ){
                U[0][i] = 1.0;
                u = -19.59745;
                p = 1000.0;
            }
            else{
                U[0][i] = 1.0;
                u = -19.59745;
                p = 0.01;
            }

            U[1][i] = U[0][i]*u;
            U[2][i] = ( p/(cts::GAMMA-1) ) + ( U[0][i]*pow(u,2)/2 );
        }

    }
    else if( cts::PROBLEM == "BLW" ){

    }
    else{
        cout << "---ERROR--- Please select correct problem ---" << endl;
    }

    return U;
}
