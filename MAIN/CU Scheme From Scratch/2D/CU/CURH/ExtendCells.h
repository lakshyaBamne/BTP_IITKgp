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

    void extend_plr(vvd& n, vvd& s, vvd& e, vvd& w, vvd& ne, vvd& nw, vvd& se, vvd& sw);
}

void ExtendCells::extend_conserved_variables(vvvd& U){
    int nx = U[0].size()-2;
    int ny = U[0][0].size()-2;
    
    if( cts::BC == "FREE" ){
        
        // extend in vertical direction
        for(int i=1 ; i<=nx ; i++){
            
            // inner loop runs 4 times for the 4 conserved variables
            for(int u=0 ; u<4 ; u++){
                U[u][i][0] = U[u][i][1];
                U[u][i][nx+1] = U[u][i][nx];    
            }

        }

        // extend in horizontal direction
        for(int j=1 ; j<=ny ; j++){
            
            // inner loop runs 4 times for the 4 conserved variables
            for(int u=0 ; u<4 ; u++){
                U[u][0][j] = U[u][1][j];
                U[u][nx+1][j] = U[u][nx][j];
            }
            
        }

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions ---" << endl;
    }
}

void ExtendCells::extend_plr(vvd& n, vvd& s, vvd& e, vvd& w, vvd& ne, vvd& nw, vvd& se, vvd& sw){

    if( cts::BC == "FREE" ){

        for(int i=1 ; i<=cts::ny ; i++){
            e[0][i] = w[1][i];
            w[cts::nx+1][i] = e[cts::nx][i];

            ne[0][i] = nw[1][i];
            nw[cts::nx+1][i] = ne[cts::nx][i];

            se[0][i] = sw[1][i];
            sw[cts::nx+1][i] = se[cts::nx][i];
        }

        for(int i=1 ; i<=cts::nx ; i++){
            n[i][0] = s[i][1];
            s[i][cts::ny+1] = n[i][cts::ny];

            sw[i][cts::ny+1] = nw[i][cts::ny];
            se[i][cts::ny+1] = ne[i][cts::ny];

            ne[i][0] = se[i][1];
            nw[i][0] = sw[i][1];
        }

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions --- " << endl;
    }

}

