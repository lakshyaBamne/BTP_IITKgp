#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<cmath>

#include "Constants.h"
#include "Utility.h"
#include "ExtendCells.h"
#include "PrimitiveVariables.h"

using namespace std;

namespace CTS = Constants;
namespace UTL = Utility;
namespace EXC = ExtendCells;
namespace PRV = PrimitiveVariables;

namespace CUNumericalFlux{
    // CU Numerical Flux (old)
    vector<vector<double>> get_cu_flux(vvd& cons_vars, pss& initial_conditions, double& dt, double& dx, double& t, pdd& time);

    // LLF Numerical Flux
    vector<vector<double>> get_cu_flux_llf(vvd& cons_vars, pss& initial_conditions, double& dt, double& dx, double& t, pdd& time);
    
    // HLL Numerical Flux
    vector<vector<double>> get_cu_flux_hll(vvd& cons_vars, pss& initial_conditions, double& dt, double& dx, double& t, pdd& time);
}

vector<vector<double>> CUNumericalFlux::get_cu_flux( vvd& cons_vars, pss& initial_conditions, double& dt, double& dx, double& t, pdd& time){
    int size = cons_vars[0].size();

    //! Step-1 Calculate the Piecewise Linear Reconstructions
    vector< vector< vector<double> > > PLR;

    for(int u=0 ; u<3 ; u++){
        double val1, val2, val3;
        double slx;

        vector< vector<double> > one_plr;
        vector<double> east(size-2);
        vector<double> west(size-2);

        for(int i=1 ; i<size-1 ; i++){
            // slope calculation
            val1 = CTS::THETA * ( cons_vars[u][i+1]-cons_vars[u][i] );
            val2 = ( cons_vars[u][i+1]-cons_vars[u][i-1] )/2;
            val3 = CTS::THETA * ( cons_vars[u][i]-cons_vars[u][i-1] );

            slx = UTL::minmod(val1, val2, val3);

            east[i-1] = cons_vars[u][i] + slx/2;
            west[i-1] = cons_vars[u][i] - slx/2;
        }

        // extend cells
        if( u == 0 ){
            EXC::extend_plr(initial_conditions.second, "density", east, west);
        }
        else if( u == 1 ){
            EXC::extend_plr(initial_conditions.second, "momentum", east, west);
        }
        else{
            EXC::extend_plr(initial_conditions.second, "energy", east, west);
        }

        one_plr.push_back(east);
        one_plr.push_back(west);

        PLR.push_back(one_plr);
    }

    /*
        * All the computations that follow are done on the N+1 Interfaces in each finite volume cell
        ! -> There are N+1 interfaces for N computational grid points 
    */
    //! Step-2 Calculate Local Speeds of Propagation
    vector<double> ap(size-1); // a+
    vector<double> am(size-1); // a-

    double uE, uW, pE, pW;
    double f1E, f2E, f3E, f1W, f2W, f3W;

    for(int i=0 ; i<size-1 ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        pE = (CTS::GAMMA-1)*(PLR[2][0][i] - 0.5*PLR[0][0][i]*pow(uE,2));
        pW = (CTS::GAMMA-1)*(PLR[2][1][i+1] - 0.5*PLR[0][1][i+1]*pow(uW,2));

        double speed_sound_E =  sqrtf( (CTS::GAMMA*pE)/PLR[0][0][i] );
        double speed_sound_W = sqrtf( (CTS::GAMMA*pW)/PLR[0][1][i+1] );

        // a+
        ap[i] = max( 0.0 , max( (uE + speed_sound_E) , (uW + speed_sound_W) ) );

        // a-
        am[i] = min( 0.0 , min( (uE - speed_sound_E) , (uW - speed_sound_W) ) );
    }

    //! Step-3 Update the Time step using CFL Conditions
    double amax=0.0;
    for(int i=0 ; i<size-1 ; i++){
        amax = max(amax, max(ap[i], -am[i]));
    }
    dt = CTS::CFL*dx/amax;

    //! Step-4 Calculate the Old CU Numerical Flux
    vector<vector<double>> cu_flux;

    vector<double> cu_f1(size);
    vector<double> cu_f2(size);
    vector<double> cu_f3(size);

    double dist, prod;

    for(int i=0 ; i<size-1 ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        pE = (CTS::GAMMA-1)*(PLR[2][0][i] - 0.5*PLR[0][0][i]*pow(uE,2));
        pW = (CTS::GAMMA-1)*(PLR[2][1][i+1] - 0.5*PLR[0][1][i+1]*pow(uW,2));
        
        f1E = PLR[1][0][i];
        f2E = PLR[1][0][i]*uE + pE;
        f3E = uE*( PLR[2][0][i] + pE );

        f1W = PLR[1][1][i+1];
        f2W = PLR[1][1][i+1]*uW + pW;
        f3W = uW*( PLR[2][1][i+1] + pW );

        dist = ap[i] - am[i];
        prod = ap[i]*am[i];

        if( dist > CTS::EPSILON ){
            // a+ - a- appears in the denominator so it cant be too small

            cu_f1[i] = ( ap[i]*f1E - am[i]*f1W + prod*( PLR[0][1][i+1] - PLR[0][0][i] ) ) / ( dist );
            cu_f2[i] = ( ap[i]*f2E - am[i]*f2W + prod*( PLR[1][1][i+1] - PLR[1][0][i] ) ) / ( dist );
            cu_f3[i] = ( ap[i]*f3E - am[i]*f3W + prod*( PLR[2][1][i+1] - PLR[2][0][i] ) ) / ( dist );
        }
        else{
            cu_f1[i] = 0.5 * ( f1E + f1W );
            cu_f2[i] = 0.5 * ( f2E + f2W );
            cu_f3[i] = 0.5 * ( f3E + f3W );
        }
    }

    cu_flux.push_back(cu_f1);
    cu_flux.push_back(cu_f2);
    cu_flux.push_back(cu_f3);

    return cu_flux;
}

vector<vector<double>> CUNumericalFlux::get_cu_flux_llf( vvd& cons_vars, pss& initial_conditions, double& dt, double& dx, double& t, pdd& time){
    int size = cons_vars[0].size();

    //! Step-1 Calculate the Piecewise Linear Reconstructions
    vector< vector< vector<double> > > PLR;

    for(int u=0 ; u<3 ; u++){
        double val1, val2, val3;
        double slx;

        vector< vector<double> > one_plr;
        vector<double> east(size-2);
        vector<double> west(size-2);

        for(int i=1 ; i<size-1 ; i++){
            val1 = cons_vars[u][i] - cons_vars[u][i-1];
            val2 = cons_vars[u][i+1] - cons_vars[u][i];

            slx = UTL::minmod(val1, val3);

            // slx = 0;

            east[i-1] = cons_vars[u][i] + slx/2;
            west[i-1] = cons_vars[u][i] - slx/2;
        }

        // extend cells
        if( u == 0 ){
            EXC::extend_plr(initial_conditions.second, "density", east, west);
        }
        else if( u == 1 ){
            EXC::extend_plr(initial_conditions.second, "momentum", east, west);
        }
        else{
            EXC::extend_plr(initial_conditions.second, "energy", east, west);
        }

        one_plr.push_back(east);
        one_plr.push_back(west);

        PLR.push_back(one_plr);
    }

    /*
        * All the computations that follow are done on the N+1 Interfaces in each finite volume cell
        ! -> There are N+1 interfaces for N computational grid points 
    */
    //! Step-2 Calculate Local Speeds of Propagation
    vector<double> ap(size-1); // a+

    vector<double> a(size-1); // alpha

    double uE, uW, pE, pW;
    double f1E, f2E, f3E, f1W, f2W, f3W;

    for(int i=0 ; i<size-1 ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        pE = (CTS::GAMMA-1)*(PLR[2][0][i] - 0.5*PLR[0][0][i]*pow(uE,2));
        pW = (CTS::GAMMA-1)*(PLR[2][1][i+1] - 0.5*PLR[0][1][i+1]*pow(uW,2));

        double speed_sound_E =  sqrtf( (CTS::GAMMA*pE)/PLR[0][0][i] );
        double speed_sound_W = sqrtf( (CTS::GAMMA*pW)/PLR[0][1][i+1] );

        // a+
        ap[i] = max( 0.0,  max( ( abs(uE) + abs(speed_sound_E) ) , ( abs(uW) + abs(speed_sound_W) ) ) );
        // ap[i] = max( ( abs(uE) + speed_sound_E)  , ( abs(uW) + speed_sound_W ) );
    
        // used to update time step
        a[i] = max( ( abs(uE) + speed_sound_E) , ( abs(uW) + speed_sound_W) );
    }

    //! Step-3 Update the Time step using CFL Conditions
    double amax=0.0;
    for(int i=0 ; i<size-1 ; i++){
        amax = max(amax, ap[i]);
        // amax = max(amax, a[i]);
    }
    dt = CTS::CFL*dx/amax;

    //! Step-4 Calculate the Old CU Numerical Flux
    vector<vector<double>> cu_flux;

    vector<double> cu_f1(size);
    vector<double> cu_f2(size);
    vector<double> cu_f3(size);

    double dist, prod;

    for(int i=0 ; i<size-1 ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        pE = (CTS::GAMMA-1)*(PLR[2][0][i] - 0.5*PLR[0][0][i]*pow(uE,2));
        pW = (CTS::GAMMA-1)*(PLR[2][1][i+1] - 0.5*PLR[0][1][i+1]*pow(uW,2));
        
        f1E = PLR[1][0][i];
        f2E = PLR[1][0][i]*uE + pE;
        f3E = uE*( PLR[2][0][i] + pE );

        f1W = PLR[1][1][i+1];
        f2W = PLR[1][1][i+1]*uW + pW;
        f3W = uW*( PLR[2][1][i+1] + pW );

        cu_f1[i] = ( ( f1E + f1W ) - ap[i]*( PLR[0][1][i+1] - PLR[0][0][i] ) ) / 2.0;
        cu_f2[i] = ( ( f2E + f2W ) - ap[i]*( PLR[1][1][i+1] - PLR[1][0][i] ) ) / 2.0;
        cu_f3[i] = ( ( f3E + f3W ) - ap[i]*( PLR[2][1][i+1] - PLR[2][0][i] ) ) / 2.0;

        // cu_f1[i] = ( ( f1E + f1W ) - a[i]*( PLR[0][1][i+1] - PLR[0][0][i] ) ) / 2.0;
        // cu_f2[i] = ( ( f2E + f2W ) - a[i]*( PLR[1][1][i+1] - PLR[1][0][i] ) ) / 2.0;
        // cu_f3[i] = ( ( f3E + f3W ) - a[i]*( PLR[2][1][i+1] - PLR[2][0][i] ) ) / 2.0;
    }

    cu_flux.push_back(cu_f1);
    cu_flux.push_back(cu_f2);
    cu_flux.push_back(cu_f3);

    return cu_flux;
}

vector<vector<double>> CUNumericalFlux::get_cu_flux_hll( vvd& cons_vars, pss& initial_conditions, double& dt, double& dx, double& t, pdd& time){
    int size = cons_vars[0].size();

    //! Step-1 Calculate the Piecewise Linear Reconstructions
    vector< vector< vector<double> > > PLR;

    for(int u=0 ; u<3 ; u++){
        double val1, val2, val3;
        double slx;

        vector< vector<double> > one_plr;
        vector<double> east(size-2);
        vector<double> west(size-2);

        for(int i=1 ; i<size-1 ; i++){
            // general minmod limiter
            val1 = CTS::THETA * ( cons_vars[u][i+1]-cons_vars[u][i] );
            val2 = ( cons_vars[u][i+1]-cons_vars[u][i-1] )/2;
            val3 = CTS::THETA * ( cons_vars[u][i]-cons_vars[u][i-1] );

            // simpler minmod limiter
            // val1 = cons_vars[u][i+1]-cons_vars[u][i];
            // val2 = cons_vars[u][i]-cons_vars[u][i-1];

            slx = UTL::minmod(val1, val2, val3);
            // slx = UTL::minmod(val1, val2);
            // slx = 0;

            east[i-1] = cons_vars[u][i] + slx/2;
            west[i-1] = cons_vars[u][i] - slx/2;
        }

        // extend cells
        if( u == 0 ){
            EXC::extend_plr(initial_conditions.second, "density", east, west);
        }
        else if( u == 1 ){
            EXC::extend_plr(initial_conditions.second, "momentum", east, west);
        }
        else{
            EXC::extend_plr(initial_conditions.second, "energy", east, west);
        }

        one_plr.push_back(east);
        one_plr.push_back(west);

        PLR.push_back(one_plr);
    }

    /*
        * All the computations that follow are done on the N+1 Interfaces in each finite volume cell
        ! -> There are N+1 interfaces for N computational grid points 
    */
    //! Step-2 Calculate Local Speeds of Propagation
    vector<double> sE(size-1);
    vector<double> sW(size-1);

    vector<double> acosW(size-1);
    vector<double> acosE(size-1);

    // double uE, uW, pE, pW;
    double f1E, f2E, f3E, f1W, f2W, f3W;

    vd uE(size);
    vd uW(size);

    vd pE(size);
    vd pW(size);

    for(int i=0 ; i<size ; i++){
        uE[i] = PLR[1][0][i] / PLR[0][0][i];
        uW[i] = PLR[1][1][i] / PLR[0][1][i];
    }

    for(int i=0 ; i<size ; i++){ 
        pE[i] = (CTS::GAMMA-1)*( PLR[2][0][i] - 0.5*PLR[0][0][i]*pow(uE[i],2) );
        pW[i] = (CTS::GAMMA-1)*( PLR[2][1][i] - 0.5*PLR[0][1][i]*pow(uW[i],2) );
    }

    for(int i=0 ; i<size-1 ; i++){
        acosE[i] = sqrtf( (CTS::GAMMA*pE[i])/PLR[0][0][i] );
        acosW[i] = sqrtf( (CTS::GAMMA*pW[i+1])/PLR[0][1][i+1] );

        sE[i] = min( uE[i] - acosE[i] , uE[i+1] - acosW[i] );
        sW[i] = min( uE[i] + acosE[i] , uW[i+1] + acosW[i] );
    }

    //! Step-3 Update the Time step using CFL Conditions
    double amax=0.0;
    for(int i=0 ; i<size-1 ; i++){
        amax = max(amax, max( abs(uE[i])+acosE[i] , abs(uW[i+1])+acosW[i] ));
    }
    dt = CTS::CFL*dx/amax;

    //! Step-4 Calculate the Old CU Numerical Flux
    vector<vector<double>> cu_flux;

    vector<double> cu_f1(size);
    vector<double> cu_f2(size);
    vector<double> cu_f3(size);

    for(int i=0 ; i<size-1 ; i++){
        f1E = PLR[1][0][i];
        f2E = PLR[1][0][i]*uE[i] + pE[i];
        f3E = uE[i]*( PLR[2][0][i] + pE[i] );

        f1W = PLR[1][1][i+1];
        f2W = PLR[1][1][i+1]*uW[i+1] + pW[i+1];
        f3W = uW[i+1]*( PLR[2][1][i+1] + pW[i+1] );

        if( sE[i]<=0 && sW[i]>=0 ){
            cu_f1[i] = sW[i]*f1E - sE[i]*f1W + sE[i]*sW[i]*( PLR[0][1][i+1] - PLR[0][0][i] );
            cu_f2[i] = sW[i]*f2E - sE[i]*f2W + sE[i]*sW[i]*( PLR[1][1][i+1] - PLR[1][0][i] );
            cu_f3[i] = sW[i]*f3E - sE[i]*f3W + sE[i]*sW[i]*( PLR[2][1][i+1] - PLR[2][0][i] );
        }
        else{
            if( sE[i]>=0 ){
                cu_f1[i] = f1E;
                cu_f2[i] = f2E;
                cu_f3[i] = f3E;
            }
            else if( sW[i]<=0 ){
                cu_f1[i] = f1W;
                cu_f2[i] = f2W;
                cu_f3[i] = f3W;
            }
        }

    }

    cu_flux.push_back(cu_f1);
    cu_flux.push_back(cu_f2);
    cu_flux.push_back(cu_f3);

    return cu_flux;
}


