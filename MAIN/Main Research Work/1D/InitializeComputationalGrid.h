/*
    ! Namespace to initialize the computational grid
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

namespace InitializeComputationalGrid{
    vector<double> make_grid();
}

vector<double> InitializeComputationalGrid::make_grid(){
    vector<double> grid(cts::n+2, 0);

    for(int i=1 ; i<=cts::n ; i++){
        grid[i] = cts::domx.first + cts::dx*( i - 0.5 );
    }

    return grid;
}


