#pragma once

#include<iostream>
#include<vector>
#include<cmath>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace PrimitiveVariables{
    vvd get_u(const vvvd& PLR);
    vvd get_p(const vvvd& PLR);
}

vvd PrimitiveVariables::get_u(const vvvd& PLR){
    vd uE(cts::n+2, 0);
    vd uW(cts::n+2, 0);

    for(int i=0 ; i<=cts::n+1 ; i++){
        uE[i] = PLR[1][0][i] / PLR[0][0][i];
        uW[i] = PLR[1][1][i] / PLR[0][1][i];
    }

    return {uE, uW};
}

vvd PrimitiveVariables::get_p(const vvvd& PLR){
    vd pE(cts::n+2, 0);
    vd pW(cts::n+2, 0);

    for(int i=0 ; i<=cts::n+1 ; i++){
        pE[i] = (cts::GAMMA-1)*( PLR[2][0][i] - ( 0.5*pow(PLR[1][0][i] , 2)/PLR[0][0][i] ) );
        pW[i] = (cts::GAMMA-1)*( PLR[2][1][i] - ( 0.5*pow(PLR[1][1][i] , 2)/PLR[0][1][i] ) );
    }

    return {pE, pW};
}

