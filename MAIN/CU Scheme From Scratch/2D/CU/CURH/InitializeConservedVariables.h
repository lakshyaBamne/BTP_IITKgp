#pragma once

#include<iostream>
#include<vector>
#include<string>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace InitializeConservedVariables{
    void initialize_conserved_variables(vvvd& U, vd gridx, vd gridy);
}

void InitializeConservedVariables::initialize_conserved_variables(vvvd& U, vd gridx, vd gridy){
    

    if( cts::PROBLEM == "MCW" ){

    }
    else if( cts::PROBLEM == "2DR" ){

        for(int i=1 ; i<=cts::nx ; i++){
            for(int j=1 ; j<=cts::ny ; j++){
                // primitive variables to be used temporarily in the loop
                double u, v, p; 

                if( gridx[i] >= 1.0 && gridy[j] >= 1.0 ){
                    U[0][i][j] = 1.5;
                    u = 0;
                    v = 0;
                    p = 1.5;
                }
                else if( gridx[i] < 1.0 && gridy[j] > 1.0 ){
                    U[0][i][j] = 0.5323;
                    u = 1.206;
                    v = 0;
                    p = 0.3;
                }
                else if( gridx[i] < 1.0 && gridy[j] < 1.0 ){
                    U[0][i][j] = 0.138;
                    u = 1.206;
                    v = 1.206;
                    p = 0.029;
                }
                else{
                    U[0][i][j] = 0.5323;
                    u = 0;
                    v = 1.206;
                    p = 0.3;
                }

                U[1][i][j] = U[0][i][j]*u; // mx
                U[2][i][j] = U[0][i][j]*v; // my
                U[3][i][j] = ( p / (cts::GAMMA-1) ) + 0.5*U[1][i][j]*u + 0.5*U[2][i][j]*v; // E

            }
        }

    }
    else if( cts::PROBLEM == "EXP" ){

    }
    else if( cts::PROBLEM == "IMP1" ){

    }
    else if( cts::PROBLEM == "IMP2" ){

    }
    else if( cts::PROBLEM == "KHI1" ){

    }
    else if( cts::PROBLEM == "KHI2" ){

    }
    else if( cts::PROBLEM == "KHI3" ){

    }
    else if( cts::PROBLEM == "RTI" ){

    }
    else{
        cout << "---ERROR--- Please enter correct problem for the simulation ---" << endl;
    }

    cout << "--- Initialized Conserved variables for => " << cts::PROBLEM << " ---" << endl;
}
