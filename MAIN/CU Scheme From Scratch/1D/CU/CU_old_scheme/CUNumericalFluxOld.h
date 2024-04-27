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
#include "InitializeConservedVariables.h"

using namespace std;

namespace cts = Constants;
namespace utl = Utility;
namespace opr = OutputResult;
namespace exc = ExtendCells;
namespace icv = InitializeConservedVariables;

namespace CUNumericalFLux{
    /*
        ! Function to calculate the CU Numerical Flux (old version) 
        ! and update time using CFL conditions
    */
    vvd get_old_cu_flux(vvd& U, bool update_time, double& dt, double T);

}

vvd CUNumericalFLux::get_old_cu_flux(vvd& U, bool update_time, double& dt, double T){
    /*
        ! First we need to calculate the Piecewise Linear Reconstructions
    */

    vd rhoE(cts::n+2, 0);
    vd rhoW(cts::n+2, 0);
    vd mxE(cts::n+2, 0);
    vd mxW(cts::n+2, 0);
    vd EE(cts::n+2, 0);
    vd EW(cts::n+2, 0);

    double slx1, slx2, slx3;
    
    for(int i=1 ; i<=cts::n ; i++){

        slx1 = utl::minmod( cts::THETA*( U[0][i] - U[0][i-1] ) , 0.5*(U[0][i+1] - U[0][i-1]) , cts::THETA*( U[0][i+1] - U[0][i] ) );
        rhoE[i] = U[0][i] + 0.5*slx1;
        rhoW[i] = U[0][i] - 0.5*slx1;

        slx2 = utl::minmod( cts::THETA*( U[1][i] - U[1][i-1] ) , 0.5*(U[1][i+1] - U[1][i-1]) , cts::THETA*( U[1][i+1] - U[1][i] ) );
        mxE[i] = U[1][i] + 0.5*slx2;
        mxW[i] = U[1][i] - 0.5*slx2;
        
        slx3 = utl::minmod( cts::THETA*( U[2][i] - U[2][i-1] ) , 0.5*(U[2][i+1] - U[2][i-1]) , cts::THETA*( U[2][i+1] - U[2][i] ) );
        EE[i] = U[2][i] + 0.5*slx3;
        EW[i] = U[2][i] - 0.5*slx3;

    }

    exc::extend_plr(rhoE, rhoW, "Density");
    exc::extend_plr(mxE, mxW, "Momentum");
    exc::extend_plr(EE, EW, "Energy");

    /*
        ! Now we can calculate the CU Numerical flux (old version)
        ! with local speeds of propagation as an intermediate

        * These are only calculate for the N-1 computational grids
    */

    vd F1(cts::n+1, 0);
    vd F2(cts::n+1, 0);
    vd F3(cts::n+1, 0);

    vd ap(cts::n+1, 0);
    // vd am(cts::n+1, 0);

    for(int i=0 ; i<=cts::n ; i++){
        // calculate primitive variables
        double uE = mxE[i] / rhoE[i];
        double uW = mxW[i+1] / rhoW[i+1];
        double pE = ( cts::GAMMA - 1 )*( EE[i] - 0.5*rhoE[i]*pow( uE,2 ) );
        double pW = ( cts::GAMMA - 1 )*( EW[i+1] - 0.5*rhoW[i+1]*pow( uW,2 ) );
    
        // local speeds of propagation
        ap[i] = max(
            0.0,
            max(
                uE + sqrt( cts::GAMMA*pE/rhoE[i] ),
                uW + sqrt( cts::GAMMA*pW/rhoW[i+1] )
            )
        );

        // am[i] = min(
        //     0.0,
        //     min(
        //         uE - sqrt( cts::GAMMA*pE/rhoE[i] ),
        //         uW - sqrt( cts::GAMMA*pW/rhoW[i+1] )
        //     )
        // );

        // flux variables
        double f1E = mxE[i];
        double f2E = rhoE[i]*pow(uE,2) + pE;
        double f3E = uE*( EE[i] + pE );

        double f1W = mxW[i+1];
        double f2W = rhoW[i+1]*pow(uW,2) + pW;
        double f3W = uW*( EW[i+1] + pW );

        double dist = ap[i] - am[i];
        double prod = ap[i]*am[i];

        if( dist > cts::EPSILON ){
            F1[i] = ( (ap[i]*f1E - am[i]*f1W) + prod*( rhoW[i+1] - rhoE[i] ) ) / dist;
            F2[i] = ( (ap[i]*f2E - am[i]*f2W) + prod*( mxW[i+1] - mxE[i] ) ) / dist;
            F3[i] = ( (ap[i]*f2E - am[i]*f2W) + prod*( EW[i+1] - EE[i] ) ) / dist;
        }
        else{
            F1[i] = 0.5*( f1E + f1W );
            F2[i] = 0.5*( f2E + f2W );
            F3[i] = 0.5*( f2E + f2W );
        }

    }

    /*
        ! Update time using CFL conditions
    */
    if( update_time ){

        double amax = 0;

        for(int i=0 ; i<=cts::n ; i++){
            amax = max( amax, max( ap[i] , -am[i] ) );
        }

        dt = cts::CFL*cts::dx/amax;

        if( T+dt > cts::time.second ){
            dt = cts::time.second - T;
        }

    }

    return {F1, F2, F3};
}
