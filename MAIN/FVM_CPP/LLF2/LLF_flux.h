#pragma once

#include "Constants.h"

using namespace std;

namespace CTS = Constants;

namespace LLF{
    // function to calculate the intermediate alpha term for LLF FLux
    vd get_alpha(const vvd& Pp, const vvd& Pm, const vvd& acos);

    // function to calculate the LLF Fluxes at the interfaces
    vvd get_llf_fluxes(const vvd& Fp, const vvd& Fm, const vd& alpha, const vvd& Up, const vvd& Um);
}

vd LLF::get_alpha(const vvd& Pp, const vvd& Pm, const vvd& acos){
    vd alpha(CTS::N+1);

    for(int i=0 ; i<=CTS::N ; i++){
        alpha[i] = max(
            abs(Pp[1][i]) + acos[0][i],
            abs(Pm[1][i]) + acos[1][i]
        );
    }

    return alpha;
}

vvd LLF::get_llf_fluxes(const vvd& Fp, const vvd& Fm, const vd& alpha, const vvd& Up, const vvd& Um){
    vvd llf_flx(3, vd(CTS::N+1));

    for(int i=0 ; i<3 ; i++){
        for(int j=0 ; j<=CTS::N ; j++){
            llf_flx[i][j] = 0.5*(Fp[i][j] + Fm[i][j]) - 0.5*alpha[j]*(Up[i][j] - Um[i][j]);
        }
    }

    return llf_flx;
}

