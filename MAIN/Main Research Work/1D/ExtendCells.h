/*
    ! Namespace contains functions to extend boundaries
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace ExtendCells{
    void extend_cells(vvd& U);
    void extend_plr(vd& east, vd& west, string var);
}

void ExtendCells::extend_cells(vvd& U){

    if( cts::BC == "FREE" ){

        for(int u=0 ; u<3 ; u++){
            U[u][0] = U[u][1];
            U[u][cts::n+1] = U[u][cts::n];
        }

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions ---" << endl;
    }

}

void ExtendCells::extend_plr(vd& east, vd& west, string var){

    if( cts::BC == "FREE" ){
        
        east[0] = west[1];
        west[cts::n+1] = east[cts::n];

    }
    else if( cts::BC == "SOLIDWALL" ){

    }
    else if( cts::BC == "PERIODIC" ){
        
    }
    else if( cts::BC == "REFLECTIVE" ){
        
    }
    else{
        cout << "---ERROR--- Please select correct boundary conditions ---" << endl;
    }

}

