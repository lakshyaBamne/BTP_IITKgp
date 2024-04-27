#pragma once

#include "Constants.h"

using namespace std;

namespace CTS = Constants;

namespace InterfaceVariables{
    // function to calculate the fluxes at the interfaces
    vvvd get_interface_fluxes(const vvd& Up, const vvd& Um);

    // function to calculate primitive variables on the interfaces
    vvvd get_primitive_variables(const vvd& Up, const vvd& Um);

    // function to calculate the acoustic speed on the interfaces
    vvd get_acoustic_speeds(const vvd& Pp, const vvd& Pm);

}

vvvd InterfaceVariables::get_interface_fluxes(const vvd& Up, const vvd& Um){
    vvd Fp(3, vd(CTS::N+1));
    vvd Fm(3, vd(CTS::N+1));

    for(int i=0 ; i<=CTS::N ; i++){
        // F(U+)
        Fp[0][i] = Up[1][i];
        Fp[1][i] = (pow(Up[1][i],2)/Up[0][i]) + (CTS::GAMMA-1)*(Up[2][i] - 0.5*(pow(Up[1][i],2)/Up[0][i]));
        Fp[2][i] = (Up[1][i]/Up[0][i])*(CTS::GAMMA*Up[2][i] - ( 0.5*(CTS::GAMMA-1)*(pow(Up[1][i],2)/Up[0][i]) ));

        // F(U-)
        Fm[0][i] = Um[1][i];
        Fm[1][i] = (pow(Um[1][i],2)/Um[0][i]) + (CTS::GAMMA-1)*(Um[2][i] - 0.5*(pow(Um[1][i],2)/Um[0][i]));
        Fm[2][i] = (Um[1][i]/Um[0][i])*(CTS::GAMMA*Um[2][i] - ( 0.5*(CTS::GAMMA-1)*(pow(Um[1][i],2)/Um[0][i]) ));
    }

    return {Fp, Fm};
}

vvvd InterfaceVariables::get_primitive_variables(const vvd& Up, const vvd& Um){
    vvd Pp(3, vd(CTS::N+1));
    vvd Pm(3, vd(CTS::N+1));

    for(int i=0 ; i<=CTS::N ; i++){
        // P+
        Pp[0][i] = Up[0][i];
        Pp[1][i] = Up[1][i]/Up[0][i];
        Pp[2][i] = (CTS::GAMMA-1)*(Up[2][i] - 0.5*(pow(Up[1][i],2)/Up[0][i]));

        // P-
        Pm[0][i] = Um[0][i];
        Pm[1][i] = Um[1][i]/Um[0][i];
        Pm[2][i] = (CTS::GAMMA-1)*(Um[2][i] - 0.5*(pow(Um[1][i],2)/Um[0][i]));
    }

    return {Pp, Pm};
}

vvd InterfaceVariables::get_acoustic_speeds(const vvd& Pp, const vvd& Pm){
    vd acos_p(CTS::N+1);
    vd acos_m(CTS::N+1);

    for(int i=0 ; i<=CTS::N ; i++){
        acos_p[i] = sqrt(CTS::GAMMA*Pp[2][i]/Pp[0][i]);
        acos_m[i] = sqrt(CTS::GAMMA*Pm[2][i]/Pm[0][i]);
    }   

    return {acos_p, acos_m}; 
}
