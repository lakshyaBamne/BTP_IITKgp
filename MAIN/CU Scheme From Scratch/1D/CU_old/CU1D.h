/*
    Functions to implement the Second Order Semi Discrete Central Upwind Scheme
    to Euler System of Gas Dynamics in 1 Dimension
*/
#pragma once

#include<iostream>
#include<vector>
#include<cmath>
#include<utility>
#include<string>

#include "ExtendCells.h"
#include "Constants.h"
#include "Parameters.h"
#include "Utility.h"

using namespace std;

namespace EXC = ExtendCells;
namespace CTS = Constants;
namespace PRM = Parameters;
namespace UTL = Utility;

#define ll long long int
#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define pdd pair<double,double>
#define pss pair<string,string>
#define pii pair<ll,ll>

namespace CentralUpwindScheme{

    // Function to calculate the CU Numerical Flux given the conserved variables
    vvd get_cu_numerical_flux(vvd cons_vars, double &dt, bool newdtneed, double dx);
    
    // Function to calculate the Piecewise Linear Reconstruction given cell averages
    vvvd get_plr(vvd cons_vars);

    void update_dt(vd ap, vd am, double &dt, double dx);
}

vvvd CentralUpwindScheme::get_plr(vvd cons_vars){
    // first we need to extend the conserved variables vector
    EXC::extend_conserved_variables(cons_vars);

    ll len = cons_vars[0].size();

    vvvd plr;
    vvd rho_plr, m_plr, e_plr;

    vd east(len-2), west(len-2);
    
    double slope;

    // rho_plr
    for(int i=1 ; i<len-1 ; i++){
        // rho_plr
        slope = UTL::minmod3(
            CTS::THETA*( cons_vars[0][i] - cons_vars[0][i-1] ),
            cons_vars[0][i+1] - cons_vars[0][i-1],
            CTS::THETA*( cons_vars[0][i+1] - cons_vars[0][i] )
        );

        east[i-1] = cons_vars[0][i] + 0.5*slope;
        west[i-1] = cons_vars[0][i] - 0.5*slope;
    }

    // extend PLR before adding
    EXC::extend_vector(east, "DENSITY");
    EXC::extend_vector(west, "DENSITY");

    rho_plr.push_back(east);
    rho_plr.push_back(west);

    // clear two elements for reuse of vector
    east.pop_back();
    east.pop_back();
    west.pop_back();
    west.pop_back();

    // m_plr
    for(int i=1 ; i<len-1 ; i++){
        // rho_plr
        slope = UTL::minmod3(
            CTS::THETA*( cons_vars[1][i] - cons_vars[1][i-1] ),
            cons_vars[1][i+1] - cons_vars[1][i-1],
            CTS::THETA*( cons_vars[1][i+1] - cons_vars[1][i] )
        );

        east[i-1] = cons_vars[1][i] + 0.5*slope;
        west[i-1] = cons_vars[1][i] - 0.5*slope;
    }

        // extend PLR before adding
    EXC::extend_vector(east, "MOMENTUM");
    EXC::extend_vector(west, "MOMENTUM");

    m_plr.push_back(east);
    m_plr.push_back(west);

    // clear two elements for reuse of vector
    east.pop_back();
    east.pop_back();
    west.pop_back();
    west.pop_back();

    // m_plr
    for(int i=1 ; i<len-1 ; i++){
        // rho_plr
        slope = UTL::minmod3(
            CTS::THETA*( cons_vars[2][i] - cons_vars[2][i-1] ),
            cons_vars[2][i+1] - cons_vars[2][i-1],
            CTS::THETA*( cons_vars[2][i+1] - cons_vars[2][i] )
        );

        east[i-1] = cons_vars[2][i] + 0.5*slope;
        west[i-1] = cons_vars[2][i] - 0.5*slope;
    }

    // extend PLR before adding
    EXC::extend_vector(east, "ENERGY");
    EXC::extend_vector(west, "ENERGY");

    e_plr.push_back(east);
    e_plr.push_back(west);

    // clear two elements for reuse of vector
    east.pop_back();
    east.pop_back();
    west.pop_back();
    west.pop_back();

    // add the vectors to the main plr
    plr.push_back(rho_plr);
    plr.push_back(m_plr);
    plr.push_back(e_plr);

    return plr;

}

vvd CentralUpwindScheme::get_cu_numerical_flux(vvd cons_vars, double &dt, bool newdtneed, double dx){
    // first get the length of the domain
    ll len = cons_vars[0].size();

    // STEP-1 : Get the Piecewise Linear Reconstruction
    vvvd plr = get_plr(cons_vars);

    vvd rho_plr = plr[0];
    vvd m_plr = plr[1];
    vvd e_plr = plr[2];

    vd ap(len+1),am(len+1);

    vvd F;
    vd F1(len+1), F2(len+1), F3(len+1);

    double uE, uW, pE, pW;
    double f1E, f1W, f2E, f2W, f3E, f3W;
    
    double ustar1, ustar2, ustar3;
    double du1, du2, du3;

    double adiff, aprod;

    for(ll i=0 ; i<=len ; i++){
        uE = m_plr[0][i]/rho_plr[0][i];
        uW = m_plr[1][i+1]/rho_plr[1][i+1];
        pE = (CTS::GAMMA-1)*(e_plr[0][i] - ( rho_plr[0][i]*pow(uE,2)/2 ));
        pW = (CTS::GAMMA-1)*(e_plr[1][i+1] - ( rho_plr[1][i+1]*pow(uW,2) ));

        f1E = m_plr[0][i];
        f1W = m_plr[1][i+1];
        f2E = rho_plr[0][i]*pow(uE,2) + pE;
        f2W = rho_plr[1][i+1]*pow(uW,2) + pW;
        f3E = uE*(e_plr[0][i] + pE);
        f3W = uW*(e_plr[1][i+1] + pW);
        
        ap[i] = max(
            0.0,
            max(
                uE + sqrtf(CTS::GAMMA*pE/rho_plr[0][i]),
                uW + sqrtf(CTS::GAMMA*pW/rho_plr[1][i+1])
            )
        );

        am[i] = min(
            0.0,
            min(
                uE - sqrtf(CTS::GAMMA*pE/rho_plr[0][i]),
                uW - sqrtf(CTS::GAMMA*pW/rho_plr[1][i+1])
            )
        );

        ustar1 = (ap[i]*rho_plr[1][i+1] - am[i]*rho_plr[0][i] - f1W + f1E) / (ap[i]-am[i]);
        ustar1 = (ap[i]*m_plr[1][i+1] - am[i]*m_plr[0][i] - f2W + f2E) / (ap[i]-am[i]);
        ustar1 = (ap[i]*e_plr[1][i+1] - am[i]*e_plr[0][i] - f3W + f3E) / (ap[i]-am[i]);

        du1 = UTL::minmod2(rho_plr[1][i+1]-ustar1 , ustar1 - rho_plr[0][i]);
        du2 = UTL::minmod2(m_plr[1][i+1]-ustar2 , ustar2 - m_plr[0][i]);
        du3 = UTL::minmod2(e_plr[1][i+1]-ustar3 , ustar3 - e_plr[0][i]);

        adiff = ap[i] - am[i];
        aprod = ap[i]*am[i];
        if( adiff > CTS::EPS ){
            F1[i] = ( ap[i]*f1E - am[i]*f1W + aprod*(rho_plr[1][i+1]-rho_plr[0][i]-du1) ) / adiff;
            F2[i] = ( ap[i]*f2E - am[i]*f2W + aprod*(m_plr[1][i+1]-m_plr[0][i]-du2) ) / adiff;
            F3[i] = ( ap[i]*f3E - am[i]*f3W + aprod*(e_plr[1][i+1]-e_plr[0][i]-du3) ) / adiff;
        }
        else{
            F1[i] = 0.5*(f1E + f1W);
            F2[i] = 0.5*(f2E + f2W);
            F3[i] = 0.5*(f3E + f3W);
        }
    }

    F.push_back(F1);
    F.push_back(F2);
    F.push_back(F3);

    if( newdtneed ){
        update_dt(ap, am, dt, dx);
    }

    return F;
}

void CentralUpwindScheme::update_dt(vd ap, vd am, double &dt, double dx){
    double amax = 0;

    ll len = ap.size();

    for(ll i=0 ; i<len ; i++){
        amax = max(
            amax,
            max(
                ap[i],
                (-1)*am[i]
            )
        );
    }

    dt = (CTS::CFL * dx) / amax;
}


