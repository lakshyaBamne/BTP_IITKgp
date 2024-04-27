/*
    ! Namespace to find the Piecewise Linear Constructions
    ! on the Finite Volume Interfaces
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>

#include "Constants.h"
#include "Utility.h"
#include "ExtendCells.h"

using namespace std;

namespace cts = Constants;
namespace utl = Utility;
namespace exc = ExtendCells;

namespace PiecewisePolynomialReconstructions{
    vvd FirstOrderPPR(const vd& U, string var_name);
    vvd SecondOrderPPR(const vd& U, string var_name);
}

vvd PiecewisePolynomialReconstructions::FirstOrderPPR(const vd& U, string var_name){
    vd east(cts::n+2, 0);
    vd west(cts::n+2, 0);

    double slx;

    for(int i=1 ; i<=cts::n ; i++){
        slx = 0;

        east[i] = U[i] + 0.5*slx;
        west[i] = U[i] - 0.5*slx;
    }

    exc::extend_plr(east, west, var_name);

    return {east,west};
}

vvd PiecewisePolynomialReconstructions::SecondOrderPPR(const vd& U, string var_name){
    vd east(cts::n+2, 0);
    vd west(cts::n+2, 0);

    double slx;

    for(int i=1 ; i<=cts::n ; i++){
        slx = utl::minmod( cts::THETA*(U[i] - U[i-1]) , 0.5*(U[i+1] - U[i-1]) , cts::THETA*( U[i+1] - U[i] ) ); 
        // slx = utl::minmod( U[i] - U[i-1] , U[i+1] - U[i] ); 

        east[i] = U[i] + 0.5*slx;
        west[i] = U[i] - 0.5*slx;
    }

    // extend PLR
    exc::extend_plr(east, west, var_name);

    return {east,west};
}

