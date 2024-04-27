#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<cmath>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace ExtendCells{
    void extend_conserved_variables(vvd& U);
    void extend_plr(vd& e, vd& w, string var);
}

void ExtendCells::extend_conserved_variables(vvd& U){
    if( cts::BC == "FREE" ){
        
        // Loop runs 3 times for the conserved variables
        for(int u=0 ; u<3 ; u++){
            U[u][0] = U[u][1];
            U[u][cts::n+1] = U[u][cts::n];
        }

    }
    else if( cts::BC == "SOLIDWALL" ){

        for(int u=0 ; u<3 ; u++){
            
            if( u!=1 ){
                U[u][0] = U[u][1];
                U[u][cts::n+1] = U[u][cts::n];
            }
            else{
                U[1][0] = -U[1][1];
                U[1][cts::n+1] = -U[1][cts::n];
            }
        }

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions ---" << endl;
    }
}

void ExtendCells::extend_plr(vd& e, vd& w, string var){
    if( cts::BC == "FREE" ){
        // e[0] = w[1];
        // w[cts::n+1] = e[cts::n];

        e[0] = e[1];
        e[cts::n+1] = e[cts::n];

        w[0] = w[1];
        w[cts::n+1] = w[cts::n];
    }
    else if( cts::BC == "SOLIDWALL" ){

        if( var != "Momentum" ){
            e[0] = w[1];
            w[cts::n+1] = e[cts::n];
        }
        else{
            e[0] = -w[1];
            w[cts::n+1] = -e[cts::n];
        }

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions ---" << endl;
    }
}
