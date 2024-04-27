#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<cmath>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace ExtendCells{
    void extend_conserved_variables(vvvd& U);

    void extend_plr(vvd& n, vvd& s, vvd& e, vvd& w, string var);
}

void ExtendCells::extend_conserved_variables(vvvd& U){
    
    if( cts::BC == "FREE" ){
        
        // extend in vertical direction
        for(int i=1 ; i<=cts::nx ; i++){
            
            // inner loop runs 4 times for the 4 conserved variables
            for(int u=0 ; u<4 ; u++){
                U[u][i][0] = U[u][i][1];
                U[u][i][cts::ny+1] = U[u][i][cts::ny];    
            }

        }

        // extend in horizontal direction
        for(int j=1 ; j<=cts::ny ; j++){
            
            // inner loop runs 4 times for the 4 conserved variables
            for(int u=0 ; u<4 ; u++){
                U[u][0][j] = U[u][1][j];
                U[u][cts::nx+1][j] = U[u][cts::nx][j];
            }
            
        }

    }
    else if( cts::BC == "REFLECTFREE" ){

        // extend in vertical direction
        for(int i=1 ; i<=cts::nx ; i++){
            
            // Density
            U[0][i][0] = U[0][i][1];
            U[0][i][cts::ny+1] = U[0][i][cts::ny];    

            // MomentumX
            U[1][i][0] = U[1][i][1];
            U[1][i][cts::ny+1] = U[1][i][cts::ny];

            // MomentumY
            U[2][i][0] = -U[2][i][1];
            U[2][i][cts::ny+1] = U[2][i][cts::ny];

            // Energy
            U[3][i][0] = U[3][i][1];
            U[3][i][cts::ny+1] = U[3][i][cts::ny];

        }

        // extend in horizontal direction
        for(int j=1 ; j<=cts::ny ; j++){
            
            // Density
            U[0][0][j] = U[0][1][j];
            U[0][cts::nx+1][j] = U[0][cts::nx][j];
            
            // MomentumX
            U[1][0][j] = -U[1][1][j];
            U[1][cts::nx+1][j] = U[1][cts::nx][j];

            // MomentumY
            U[2][0][j] = U[2][1][j];
            U[2][cts::nx+1][j] = U[2][cts::nx][j];

            // Energy
            U[3][0][j] = U[3][1][j];
            U[3][cts::nx+1][j] = U[3][cts::nx][j];

        }

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions ---" << endl;
    }
}

void ExtendCells::extend_plr(vvd& n, vvd& s, vvd& e, vvd& w, string var){

    if( cts::BC == "FREE" ){

        for(int i=1 ; i<=cts::ny ; i++){
            e[0][i] = w[1][i];
            w[cts::nx+1][i] = e[cts::nx][i];
        }

        for(int i=1 ; i<=cts::nx ; i++){
            n[i][0] = s[i][1];
            s[i][cts::ny+1] = n[i][cts::ny];
        }

    }
    else if( cts::BC == "REFLECTFREE" ){

        for(int i=1 ; i<=cts::ny ; i++){

            if( var != "MomentumX" ){
                e[0][i] = w[1][i];
                w[cts::nx+1][i] = e[cts::nx][i];
            }
            else{
                e[0][i] = -w[1][i];
                w[cts::nx+1][i] = e[cts::nx][i];
            }

        }

        for(int i=1 ; i<=cts::nx ; i++){

            if( var != "MomentumY" ){
                n[i][0] = s[i][1];
                s[i][cts::ny+1] = n[i][cts::ny];
            }
            else{
                n[i][0] = -s[i][1];
                s[i][cts::ny+1] = n[i][cts::ny];
            }

        }

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions --- " << endl;
    }

}

