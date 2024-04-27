#pragma once

#include "Constants.h"

using namespace std;

namespace CTS = Constants;

namespace TimeStepCalculation{
    // function to calculate the next time step using CFL conditions
    void update_time(double& dt, const vd& u, const vd& a);
}

void TimeStepCalculation::update_time(double& dt, const vd& u, const vd& a){
    double amax;

    for(int i=1 ; i<=CTS::N ; i++){
        amax = max(
            amax,
            abs(u[i]) + a[i]
        );
    }

    dt = CTS::CFL*CTS::dx/amax;
}
