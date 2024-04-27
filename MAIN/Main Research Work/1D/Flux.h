/*
    ! Namespace contains the functions used to calculate
    ! the CentralUpwind Numerical Flux
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>

#include "Constants.h"
#include "PrimitiveVariables.h"
#include "PiecewisePolynomialReconstructions.h"

using namespace std;

namespace cts = Constants;
namespace prv = PrimitiveVariables;
namespace ppr = PiecewisePolynomialReconstructions;

namespace Flux{
    // Local Lax Friedrichs Flux
    vvd cu_flux_llf(const vvd& U, bool update_time, double& dt, double T);

    // variants of the Central Upwind Flux
    vvd cu_flux_old(const vvd& U, bool update_time, double& dt, double T);
    vvd cu_flux(const vvd& U, bool update_time, double& dt, double T);
    vvd curh_flux(const vvd& U, bool update_time, double& dt, double T);

    // HLL FLux
    vvd cu_flux_hll(const vvd& U, bool update_time, double& dt, double T);
}

vvd Flux::cu_flux_llf(const vvd& U, bool update_time, double& dt, double T){

    vvvd PLR;

    //! First we calculate the approximations
    if( cts::APPROX == "FIRST" ){
        for(int u=0 ; u<3 ; u++){
            PLR.push_back( ppr::FirstOrderPPR(U[0], "Density") );
            PLR.push_back( ppr::FirstOrderPPR(U[1], "Momentum") );
            PLR.push_back( ppr::FirstOrderPPR(U[2], "Energy") );
        }
    }
    else if( cts::APPROX == "SECOND" ){
        for(int u=0 ; u<3 ; u++){
            PLR.push_back( ppr::SecondOrderPPR(U[0], "Density") );
            PLR.push_back( ppr::SecondOrderPPR(U[1], "Momentum") );
            PLR.push_back( ppr::SecondOrderPPR(U[2], "Energy") );
        }
    }
    else{
        cout << "---ERROR--- Please enter the correct approximation type ---" << endl;
    }

    //! Now we can calculate the LLF Numerical flux
    double ap, am;
    double uE, uW, pE, pW;
    double f1E, f2E, f3E, f1W, f2W, f3W;

    vector<double> alpha(cts::n+2, 0);

    vector< vector<double> > F(3, vector<double>(cts::n+2,0));

    vvd u = prv::get_u(PLR);
    vvd p = prv::get_p(PLR);

    for(int i=0 ; i<=cts::n ; i++){
        // uE = PLR[1][0][i] / PLR[0][0][i];
        // uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        // pE = (cts::GAMMA-1)*( PLR[2][0][i] - PLR[0][0][i]*pow(uE,2)*0.5 );
        // pW = (cts::GAMMA-1)*( PLR[2][1][i+1] - PLR[0][1][i+1]*pow(uW,2)*0.5 );

        ap = sqrtf( cts::GAMMA*pW / PLR[0][1][i+1] );
        am = sqrtf( cts::GAMMA*pE / PLR[0][0][i] );

        alpha[i] = max( abs(uE)+am , abs(uW)+ap );

        f1E = PLR[1][0][i];
        // f2E = PLR[1][0][i]*uE + pE;
        // f3E = uE*(PLR[2][0][i] + pE);

        f1E = 

        f1W = PLR[1][1][i+1];
        f2W = PLR[1][1][i+1]*uW + pW;
        f3W = uW*(PLR[2][1][i+1] + pW);
    
        F[0][i] = ( f1E + f1W - alpha[i]*( PLR[0][0][i] - PLR[0][1][i+1] ) ) / 2.0;
        F[1][i] = ( f2E + f2W - alpha[i]*( PLR[1][0][i] - PLR[1][1][i+1] ) ) / 2.0;
        F[2][i] = ( f3E + f3W - alpha[i]*( PLR[2][0][i] - PLR[2][1][i+1] ) ) / 2.0;
    }

    if( update_time ){
        double amax = 0.0;

        for(int i=0 ; i<=cts::n ; i++){
            amax = max( amax , alpha[i] );
        }

        dt = cts::CFL*cts::dx / amax;
    }

    return F;
}


vvd Flux::cu_flux_old(const vvd& U, bool update_time, double& dt, double T){

    vvvd PLR;

    //! First we calculate the approximations
    if( cts::APPROX == "FIRST" ){
        for(int u=0 ; u<U.size() ; u++){
            PLR.push_back( ppr::FirstOrderPPR(U[0], "Density") );
            PLR.push_back( ppr::FirstOrderPPR(U[1], "Momentum") );
            PLR.push_back( ppr::FirstOrderPPR(U[2], "Energy") );
        }
    }
    else if( cts::APPROX == "SECOND" ){
        for(int u=0 ; u<U.size() ; u++){
            PLR.push_back( ppr::SecondOrderPPR(U[0], "Density") );
            PLR.push_back( ppr::SecondOrderPPR(U[1], "Momentum") );
            PLR.push_back( ppr::SecondOrderPPR(U[2], "Energy") );
        }
    }
    else{
        cout << "---ERROR--- Please enter the correct approximation type ---" << endl;
    }

}

vvd Flux::cu_flux(const vvd& U, bool update_time, double& dt, double T){

}

vvd Flux::curh_flux(const vvd& U, bool update_time, double& dt, double T){

}

vvd Flux::cu_flux_hll(const vvd& U, bool update_time, double& dt, double T){

}

