#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>
#include<sstream>
#include<fstream>

#include "Constants.h"
#include "Utility.h"
#include "OutputResult.h"
#include "ExtendCells.h"
#include "CUNumericalFluxOld.h"
#include "InitializeConservedVariables.h"

using namespace std;

namespace cts = Constants;
namespace utl = Utility;
namespace opr = OutputResult;
namespace exc = ExtendCells;
namespace icv = InitializeConservedVariables;
namespace cuf = CUNumericalFLux;

namespace SSPRK{
    /*
        ! Function to update the conserved variables for the next iteration
        ! using Three stage SSPRK Time Integration Scheme
    */
    void update_conserved_variables(vvd& U, double& dt, double T);
}

void SSPRK::update_conserved_variables(vvd& U, double& dt, double T){
    //! Stage-1 (if only stage-1 is used we have the Euler Forward difference integration)
    vvd F = cuf::get_old_cu_flux(U, true, dt, T);

    // calculate LAMBDA after updating time step used in SSPRK stages
    double LAMBDA = dt / cts::dx;

    vvd U1(3, vd(cts::n+2, 0));

    for(int u=0 ; u<3 ; u++){
        for(int i=1 ; i<=cts::n ; i++){
            U1[u][i] = U[u][i] - LAMBDA*( F[u][i] - F[u][i-1] );
        }
    }

    exc::extend_conserved_variables(U1);

    //! Stage-2
    vvd F1 = cuf::get_old_cu_flux(U1, false, dt, T);

    vvd U2(3, vd(cts::n+2, 0));

    for(int u=0 ; u<3 ; u++){
        for(int i=1 ; i<=cts::n ; i++){
            U2[u][i] = ( 3.0*U[u][i] + U1[u][i] - LAMBDA*( F1[u][i] - F1[u][i-1] ) ) / 4.0;
        }
    }

    exc::extend_conserved_variables(U2);

    //! Stage-3
    vvd F2 = cuf::get_old_cu_flux(U2, false, dt, T);

    vvd UN(3, vd(cts::n+2, 0));

    for(int u=0 ; u<3 ; u++){
        for(int i=1 ; i<=cts::n ; i++){
            UN[u][i] = ( U[u][i] + 2.0*U2[u][i] - 2.0*LAMBDA*( F2[u][i] - F2[u][i-1] ) ) / 3.0;
        }
    }

    exc::extend_conserved_variables(UN);

    // update conserved variables after three stages
    for(int u=0 ; u<3 ; u++){
        for(int i=0 ; i<=cts::n+1 ; i++){
            U[u][i] = UN[u][i];
        }
    }
}


