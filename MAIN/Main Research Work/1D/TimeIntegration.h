/*
    ! Namespace contains parameters used in the program for simulations
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>

#include "Constants.h"
#include "Flux.h"

using namespace std;

namespace cts = Constants;
namespace flx = Flux;

/*
    ! 3-stage Strong stability preserving Runge-Kutta time integration
*/
namespace SSPRK{
    void update_conserved_variables(vvd& U, double& dt, double T);
}

void SSPRK::update_conserved_variables(vvd& U, double& dt, double T){

    //! Stage-1
    vvd F1;
    
    if( cts::FLUX == "CU-OLD" ){
        F1 = flx::cu_flux_old(U, true, dt, T);
    }
    else if( cts::FLUX == "CU" ){
        F1 = flx::cu_flux(U, true, dt, T);
    }
    else if( cts::FLUX == "CURH" ){
        F1 = flx::curh_flux(U, true, dt, T);
    }
    else if( cts::FLUX == "LLF" ){
        F1 = flx::cu_flux_llf(U, true, dt, T);
    }
    else if( cts::FLUX == "HLL" ){
        F1 = flx::cu_flux_hll(U, true, dt, T);
    }
    else{
        cout << "---ERROR--- Please select a valid flux to be used ---" << endl;
    }

    vvd U1(U.size(), vd(cts::n+2, 0));

    for(int u=0 ; u<U.size() ; u++){
        for(int i=1 ; i<=cts::n ; i++){
            U1[u][i] = U[u][i] - (dt/cts::dx)*( F1[u][i] - F1[u][i-1] );
        }
    }

    exc::extend_cells(U1);

    //! Stage-2
    vvd F2;

    if( cts::FLUX == "CU-OLD" ){
        F2 = flx::cu_flux_old(U1, false, dt, T);
    }
    else if( cts::FLUX == "CU" ){
        F2 = flx::cu_flux(U1, false, dt, T);
    }
    else if( cts::FLUX == "CURH" ){
        F2 = flx::curh_flux(U1, false, dt, T);
    }
    else if( cts::FLUX == "LLF" ){
        F2 = flx::cu_flux_llf(U1, false, dt, T);
    }
    else if( cts::FLUX == "HLL" ){
        F2 = flx::cu_flux_hll(U1, false, dt, T);
    }
    else{
        cout << "---ERROR--- Please select a valid flux to be used ---" << endl;
    }

    vvd U2(U.size(), vd(cts::n+2, 0));

    for(int u=0 ; u<U.size() ; u++){
        for(int i=1 ; i<=cts::n ; i++){
            U2[u][i] = ( 3*U[u][i] + U1[u][i] - (dt/cts::dx)*( F2[u][i] - F2[u][i-1] ) ) / 4;
        }
    }

    exc::extend_cells(U2);

    //! Stage-3
    vvd F3;

    if( cts::FLUX == "CU-OLD" ){
        F3 = flx::cu_flux_old(U2, false, dt, T);
    }
    else if( cts::FLUX == "CU" ){
        F3 = flx::cu_flux(U2, false, dt, T);
    }
    else if( cts::FLUX == "CURH" ){
        F3 = flx::curh_flux(U2, false, dt, T);
    }
    else if( cts::FLUX == "LLF" ){
        F3 = flx::cu_flux_llf(U2, false, dt, T);
    }
    else if( cts::FLUX == "HLL" ){
        F3 = flx::cu_flux_hll(U2, false, dt, T);
    }
    else{
        cout << "---ERROR--- Please select a valid flux to be used ---" << endl;
    }

    vvd UN(U.size(), vd(cts::n+2, 0));

    for(int u=0 ; u<U.size() ; u++){
        for(int i=1 ; i<=cts::n ; i++){
            UN[u][i] = ( U[u][i] + 2*U2[u][i] - 2*(dt/cts::dx)*( F3[u][i] - F3[u][i-1] ) ) / 3;
        }
    }

    exc::extend_cells(UN);

    for(int u=0 ; u<U.size() ; u++){
        for(int i=0 ; i<=cts::n+1 ; i++){
            U[u][i] = UN[u][i];
        }
    }

}

/*
    ! Euler forward difference time integration
*/
namespace EFD{
    void update_conserved_variables(vvd& U, double& dt, double T);
}

void EFD::update_conserved_variables(vvd& U, double& dt, double T){
    //! select the correct flux to be used to simulate
    vvd F;
    
    if( cts::FLUX == "CU-OLD" ){
        F = flx::cu_flux_old(U, true, dt, T);
    }
    else if( cts::FLUX == "CU" ){
        F = flx::cu_flux(U, true, dt, T);
    }
    else if( cts::FLUX == "CURH" ){
        F = flx::curh_flux(U, true, dt, T);
    }
    else if( cts::FLUX == "LLF" ){
        F = flx::cu_flux_llf(U, true, dt, T);
    }
    else if( cts::FLUX == "HLL" ){
        F = flx::cu_flux_hll(U, true, dt, T);
    }
    else{
        cout << "---ERROR--- Please select a valid flux to be used ---" << endl;
    }

    for(int u=0 ; u<U.size() ; u++){
        for(int i=1 ; i<=cts::n ; i++){
            U[u][i] = U[u][i] - (dt/cts::dx)*( F[u][i] - F[u][i-1] );
        }
    }

    // extend cells
    exc::extend_cells(U);
}   
