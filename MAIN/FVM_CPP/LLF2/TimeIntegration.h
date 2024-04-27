#pragma once

#include "Constants.h"

using namespace std;

namespace CTS = Constants;

namespace TimeIntegrationSchemes{
    // function to implement Euler Forward Difference time integration
    void EFD(vvd& U, const vvd& llf, double dt);

    // function to update other variables
    void update_variables(const vvd& U, vvd& P, vd& a);
}

void TimeIntegrationSchemes::EFD(vvd& U, const vvd& llf, double dt){
    for(int i=0 ; i<3 ; i++){
        for(int j=1 ; j<=CTS::N ; j++){
            U[i][j] = U[i][j] - (dt/CTS::dx)*(llf[i][j] - llf[i][j-1]);
        }
    }
}

void TimeIntegrationSchemes::update_variables(const vvd& U, vvd& P, vd& a){
    for(int i=1 ; i<=CTS::N ; i++){
        P[0][i] = U[0][i];
        P[1][i] = U[1][i]/U[0][i];
        P[2][i] = (CTS::GAMMA-1)*(U[2][i] - 0.5*(pow(U[1][i],2)/U[0][i]));
    
        a[i] = sqrt(CTS::GAMMA*P[2][i]/P[0][i]);
    }
}
