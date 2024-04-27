#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<cmath>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace InitializeConservedVariables{
    void initialize_conserved_variables(vvvd& U, vd gridx, vd gridy);
}

void InitializeConservedVariables::initialize_conserved_variables(vvvd& U, vd gridx, vd gridy){
    if( cts::PROBLEM == "MCW" ){

        for(int i=1 ; i<=cts::nx ; i++){
            for(int j=1 ; j<=cts::ny ; j++){
                
                bool c1 = gridx[i]<0.1 && gridx[i]>-0.1 && gridy[j]<0.02 && gridy[j]>0;
                bool c2 = gridx[i]<0.02 && gridx[i]>-0.02 && gridy[j]>0.02 && gridy[j]<0.1;
                bool c3 = pow( gridx[i]+0.02 , 2 ) + pow( gridy[j]-0.02 , 2 ) <= pow( 0.08 , 2 );
                bool c4 = pow( gridx[i]-0.02 , 2 ) + pow( gridy[j]-0.02 , 2 ) <= pow( 0.08 , 2 );

                if( c1 || c2 || c3 || c4 ){
                    U[0][i][j] = 1.4;
                    U[1][i][j] = 0;
                    U[2][i][j] = 0.28;
                    U[3][i][j] = ( 1.0/(cts::GAMMA-1) ) + 0.5*( pow(U[1][i][j],2) + pow(U[2][i][j],2) )/U[0][i][j];
                }
                else{
                    U[0][i][j] = 1.0;
                    U[1][i][j] = 0;
                    U[2][i][j] = 0.2;
                    U[3][i][j] = ( 1.0/(cts::GAMMA-1) ) + 0.5*( pow(U[1][i][j],2) + pow(U[2][i][j],2) )/U[0][i][j];
                }
            
            }
        }

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

        for(int i=1 ; i<=cts::nx ; i++){
            for(int j=1 ; j<=cts::ny ; j++){
                
                if( pow(gridx[i],2)+pow(gridy[j],2) < 0.16 ){
                    U[0][i][j] = 1.0;
                    U[1][i][j] = 0.0;
                    U[2][i][j] = 0.0;
                    U[3][i][j] = 1.0/(cts::GAMMA-1);
                }
                else{
                    U[0][i][j] = 0.125;
                    U[1][i][j] = 0.0;
                    U[2][i][j] = 0.0;
                    U[3][i][j] = 0.1/(cts::GAMMA-1);
                }
        
            }
        }

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
