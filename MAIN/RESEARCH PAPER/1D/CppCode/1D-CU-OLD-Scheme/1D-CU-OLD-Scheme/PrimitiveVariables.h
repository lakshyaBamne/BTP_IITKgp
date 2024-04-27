#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<cmath>

#include "Constants.h"
#include "Utility.h"
#include "ExtendCells.h"

using namespace std;

namespace CTS = Constants;
namespace UTL = Utility;
namespace EXC = ExtendCells;

namespace PrimitiveVariables{
    vvd get_u(const vvvd& PLR);
    vvd get_p(const vvvd& PLR);
}

vvd PrimitiveVariables::get_u(const vvvd& PLR){
    int size = PLR[0][0].size();
    
    vd uE(size-1);
    vd uW(size-1);

    for(int i=0 ; i<size ; i++){
        uE[i] = PLR[1][0][i] / PLR[0][0][i];
        uW[i] = PLR[1][1][i] / PLR[0][1][i];
    }

    return {uE, uW};
}

vvd PrimitiveVariables::get_p(const vvvd& PLR){
    int size = PLR[0][0].size();
    
    vd pE(size-1);
    vd pW(size-1);

    double uE, uW;

    for(int i=0 ; i<size ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i] / PLR[0][1][i];
        
        pE[i] = (CTS::GAMMA-1)*( PLR[2][0][i] - 0.5*PLR[0][0][i]*pow(uE,2) );
        pW[i] = (CTS::GAMMA-1)*( PLR[2][1][i] - 0.5*PLR[0][1][i]*pow(uW,2) );
    }

    return {pE, pW};
}


