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

namespace CUNumericalFlux{
    //! 1D code
    vector< vector<double> > get_cu_flux(vector< vector<double> >& cons_vars, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time);

    //! 2D code
    vector< vector<double> > get_cu_flux_horizontal(vector< vector<double> >& cons_vars, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time, double& amax_one);
    vector< vector<double> > get_cu_flux_vertical(vector< vector<double> >& cons_vars, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time, double& amax_one);
}

//! 2D code
vector< vector<double> > CUNumericalFlux::get_cu_flux_horizontal(vector< vector<double> >& cons_vars, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time, double& amax_one){
    int size = cons_vars[0].size();

    //! Step-1 Calculate the Piecewise Linear Reconstructions
    vector< vector< vector<double> > > PLR;

    for(int u=0 ; u<4 ; u++){
        double val1, val2, val3;
        double slx;

        vector< vector<double> > one_plr;
        vector<double> east(size-2);
        vector<double> west(size-2);

        for(int i=1 ; i<size-1 ; i++){
            // slope calculation
            double val1 = CTS::THETA * ( cons_vars[u][i+1]-cons_vars[u][i] );
            double val2 = ( cons_vars[u][i+1]-cons_vars[u][i-1] )/2;
            double val3 = CTS::THETA * ( cons_vars[u][i]-cons_vars[u][i-1] );

            slx = UTL::minmod(val1, val2, val3);

            east[i-1] = cons_vars[u][i] + slx/2;
            west[i-1] = cons_vars[u][i] - slx/2;
        }

        // extend cells
        if( u == 0 ){
            EXC::extend_plr_horizontal(initial_conditions.second, "density", east, west);
        }
        else if( u == 1 ){
            EXC::extend_plr_horizontal(initial_conditions.second, "momentumX", east, west);
        }
        else if( u==2 ){
            EXC::extend_plr_horizontal(initial_conditions.second, "momentumY", east, west);
        }
        else{
            EXC::extend_plr_horizontal(initial_conditions.second, "energy", east, west);
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

    double uE, uW, vE, vW, pE, pW;
    double f1E, f2E, f3E, f4E, f1W, f2W, f3W, f4W;

    for(int i=0 ; i<size-1 ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        vE = PLR[2][0][i] / PLR[0][0][i];
        vW = PLR[2][1][i+1] / PLR[0][1][i+1];
        pE = (CTS::GAMMA-1)*(PLR[3][0][i] - 0.5*PLR[0][0][i]*pow(uE,2) - 0.5*PLR[0][0][i]*pow(vE,2));
        pW = (CTS::GAMMA-1)*(PLR[3][1][i+1] - 0.5*PLR[0][1][i+1]*pow(uW,2) - 0.5*PLR[0][1][i+1]*pow(vW,2));

        double speed_sound_E =  sqrtf( (CTS::GAMMA*pE)/PLR[0][0][i] );
        double speed_sound_W = sqrtf( (CTS::GAMMA*pW)/PLR[0][1][i+1] );

        // a+
        ap[i] = max( 0.0 , max( (uE + speed_sound_E) , (uW + speed_sound_W) ) );

        // a-
        am[i] = min( 0.0 , min( (uE - speed_sound_E) , (uW - speed_sound_W) ) );
    }

    //! Step-3 Update the Time step using CFL Conditions
    for(int i=0 ; i<size-1 ; i++){
        amax_one = max(amax_one, max(ap[i], -am[i]));
    }

    //! Step-4 Calculate the Old CU Numerical Flux
    vector<vector<double>> cu_flux;

    vector<double> cu_f1(size-1);
    vector<double> cu_f2(size-1);
    vector<double> cu_f3(size-1);
    vector<double> cu_f4(size-1);
    
    double dist, prod;

    for(int i=0 ; i<size-1 ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        vE = PLR[2][0][i] / PLR[0][0][i];
        vW = PLR[2][1][i+1] / PLR[0][1][i+1];
        pE = (CTS::GAMMA-1)*(PLR[3][0][i] - 0.5*PLR[0][0][i]*pow(uE,2) - 0.5*PLR[0][0][i]*pow(vE,2));
        pW = (CTS::GAMMA-1)*(PLR[3][1][i+1] - 0.5*PLR[0][1][i+1]*pow(uW,2) - 0.5*PLR[0][1][i+1]*pow(vW,2));
        
        f1E = PLR[1][0][i]; // rho*u
        f2E = PLR[1][0][i]*uE + pE; // rho*u*u + p
        f3E = PLR[1][0][i]*vE; // rho*u*v
        f4E = uE*( PLR[3][0][i] + pE ); // u*(E + p)

        f1W = PLR[1][1][i+1]; // rho*u
        f2W = PLR[1][1][i+1]*uW + pW; // rho*u*u + p
        f3W = PLR[1][1][i+1]*vW; // rho*u*v
        f4W = uW*( PLR[3][1][i+1] + pW ); // u*(E + p)

        dist = ap[i] - am[i];
        prod = ap[i]*am[i];

        if( dist > CTS::EPSILON ){
            // a+ - a- appears in the denominator so it cant be too small

            cu_f1[i] = ( ap[i]*f1E - am[i]*f1W + prod*( PLR[0][1][i+1] - PLR[0][0][i] ) ) / ( dist );
            cu_f2[i] = ( ap[i]*f2E - am[i]*f2W + prod*( PLR[1][1][i+1] - PLR[1][0][i] ) ) / ( dist );
            cu_f3[i] = ( ap[i]*f3E - am[i]*f3W + prod*( PLR[2][1][i+1] - PLR[2][0][i] ) ) / ( dist );
            cu_f4[i] = ( ap[i]*f4E - am[i]*f4W + prod*( PLR[3][1][i+1] - PLR[3][0][i] ) ) / ( dist );
        }
        else{
            cu_f1[i] = 0.5 * ( f1E + f1W );
            cu_f2[i] = 0.5 * ( f2E + f2W );
            cu_f3[i] = 0.5 * ( f3E + f3W );
            cu_f4[i] = 0.5 * ( f4E + f4W );
        }
    }

    cu_flux.push_back(cu_f1);
    cu_flux.push_back(cu_f2);
    cu_flux.push_back(cu_f3);
    cu_flux.push_back(cu_f4);

    return cu_flux;
}

vector< vector<double> > CUNumericalFlux::get_cu_flux_vertical(vector< vector<double> >& cons_vars, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time, double& bmax_one){
    int size = cons_vars[0].size();

    //! Step-1 Calculate the Piecewise Linear Reconstructions
    vector< vector< vector<double> > > PLR;

    for(int u=0 ; u<4 ; u++){
        double val1, val2, val3;
        double slx;

        vector< vector<double> > one_plr;
        vector<double> east(size-2);
        vector<double> west(size-2);

        for(int i=1 ; i<size-1 ; i++){
            // slope calculation
            double val1 = CTS::THETA * ( cons_vars[u][i+1]-cons_vars[u][i] );
            double val2 = ( cons_vars[u][i+1]-cons_vars[u][i-1] )/2;
            double val3 = CTS::THETA * ( cons_vars[u][i]-cons_vars[u][i-1] );

            slx = UTL::minmod(val1, val2, val3);

            east[i-1] = cons_vars[u][i] + slx/2;
            west[i-1] = cons_vars[u][i] - slx/2;
        }

        // extend cells
        if( u == 0 ){
            EXC::extend_plr_vertical(initial_conditions.second, "density", east, west);
        }
        else if( u == 1 ){
            EXC::extend_plr_vertical(initial_conditions.second, "momentumX", east, west);
        }
        else if( u==2 ){
            EXC::extend_plr_vertical(initial_conditions.second, "momentumY", east, west);
        }
        else{
            EXC::extend_plr_vertical(initial_conditions.second, "energy", east, west);
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

    double uE, uW, vE, vW, pE, pW;
    double f1E, f2E, f3E, f4E, f1W, f2W, f3W, f4W;

    for(int i=0 ; i<size-1 ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        vE = PLR[2][0][i] / PLR[0][0][i];
        vW = PLR[2][1][i+1] / PLR[0][1][i+1];
        pE = (CTS::GAMMA-1)*(PLR[3][0][i] - 0.5*PLR[0][0][i]*pow(uE,2) - 0.5*PLR[0][0][i]*pow(vE,2));
        pW = (CTS::GAMMA-1)*(PLR[3][1][i+1] - 0.5*PLR[0][1][i+1]*pow(uW,2) - 0.5*PLR[0][1][i+1]*pow(vW,2));

        double speed_sound_E =  sqrtf( (CTS::GAMMA*pE)/PLR[0][0][i] );
        double speed_sound_W = sqrtf( (CTS::GAMMA*pW)/PLR[0][1][i+1] );

        // a+
        ap[i] = max( 0.0 , max( (uE + speed_sound_E) , (uW + speed_sound_W) ) );

        // a-
        am[i] = min( 0.0 , min( (uE - speed_sound_E) , (uW - speed_sound_W) ) );
    }

    //! Step-3 Update the Time step using CFL Conditions
    for(int i=0 ; i<size-1 ; i++){
        bmax_one = max(bmax_one, max(ap[i], -am[i]));
    }

    //! Step-4 Calculate the Old CU Numerical Flux
    vector<vector<double>> cu_flux;

    vector<double> cu_f1(size-1);
    vector<double> cu_f2(size-1);
    vector<double> cu_f3(size-1);
    vector<double> cu_f4(size-1);
    
    double dist, prod;

    for(int i=0 ; i<size-1 ; i++){
        uE = PLR[1][0][i] / PLR[0][0][i];
        uW = PLR[1][1][i+1] / PLR[0][1][i+1];
        vE = PLR[2][0][i] / PLR[0][0][i];
        vW = PLR[2][1][i+1] / PLR[0][1][i+1];
        pE = (CTS::GAMMA-1)*(PLR[3][0][i] - 0.5*PLR[0][0][i]*pow(uE,2) - 0.5*PLR[0][0][i]*pow(vE,2));
        pW = (CTS::GAMMA-1)*(PLR[3][1][i+1] - 0.5*PLR[0][1][i+1]*pow(uW,2) - 0.5*PLR[0][1][i+1]*pow(vW,2));
        
        f1E = PLR[2][0][i]; // rho*v
        f2E = PLR[2][0][i]*uE; // rho*u*v
        f3E = PLR[2][0][i]*vE + pE; // rho*v*v + p
        f4E = vE*( PLR[3][0][i] + pE ); // v*(E + p)

        f1W = PLR[2][1][i+1]; // rho*v
        f2W = PLR[2][1][i+1]*uW; // rho*u*v
        f3W = PLR[2][1][i+1]*vW + pW; // rho*v*v + p
        f4W = vW*( PLR[3][1][i+1] + pW ); // v*(E + p)

        dist = ap[i] - am[i];
        prod = ap[i] * am[i];

        if( dist > CTS::EPSILON ){
            // a+ - a- appears in the denominator so it cant be too small

            cu_f1[i] = ( ap[i]*f1E - am[i]*f1W + prod*( PLR[0][1][i+1] - PLR[0][0][i] ) ) / ( dist );
            cu_f2[i] = ( ap[i]*f2E - am[i]*f2W + prod*( PLR[1][1][i+1] - PLR[1][0][i] ) ) / ( dist );
            cu_f3[i] = ( ap[i]*f3E - am[i]*f3W + prod*( PLR[2][1][i+1] - PLR[2][0][i] ) ) / ( dist );
            cu_f4[i] = ( ap[i]*f4E - am[i]*f4W + prod*( PLR[3][1][i+1] - PLR[3][0][i] ) ) / ( dist );
        }
        else{
            cu_f1[i] = 0.5 * ( f1E + f1W );
            cu_f2[i] = 0.5 * ( f2E + f2W );
            cu_f3[i] = 0.5 * ( f3E + f3W );
            cu_f4[i] = 0.5 * ( f4E + f4W );
        }
    }

    cu_flux.push_back(cu_f1);
    cu_flux.push_back(cu_f2);
    cu_flux.push_back(cu_f3);
    cu_flux.push_back(cu_f4);

    return cu_flux;
}


//! 1D code
vector< vector<double> > CUNumericalFlux::get_cu_flux( vector< vector<double> >& cons_vars, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time){
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
            double val1 = CTS::THETA * ( cons_vars[u][i+1]-cons_vars[u][i] );
            double val2 = ( cons_vars[u][i+1]-cons_vars[u][i-1] )/2;
            double val3 = CTS::THETA * ( cons_vars[u][i]-cons_vars[u][i-1] );

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

    vector<double> cu_f1(size-1);
    vector<double> cu_f2(size-1);
    vector<double> cu_f3(size-1);
    
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
