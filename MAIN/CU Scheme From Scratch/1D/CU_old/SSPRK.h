#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<utility>

#include "Constants.h"
#include "Parameters.h"
#include "CU1D.h"

using namespace std;

namespace CTS = Constants;
namespace PRM = Parameters;
namespace CUS = CentralUpwindScheme;

#define ll long long int
#define vd vector<double>
#define vvd vector< vector<double> >
#define pdd pair<double,double>
#define pss pair<string,string>
#define pii pair<ll,ll>

namespace SSPRK{
    // Function to update conserved variables for the next time step
    // given the CU Numerical Flux
    void update_conserved_variables(vvd &cons_vars, double &dt, double dx);
}

void SSPRK::update_conserved_variables(vvd &cons_vars, double &dt, double dx){

    // length of the conserved variables vector
    ll len = cons_vars[0].size();

    vvd F1 = CUS::get_cu_numerical_flux(cons_vars, dt, true, dx);

    vvd u1;
    vd u11(len), u12(len), u13(len);

    for(ll i=0 ; i<len ; i++){
        u11[i] = cons_vars[0][i] - (dt/dx)*(F1[0][i+1] - F1[0][i]);
        u12[i] = cons_vars[1][i] - (dt/dx)*(F1[1][i+1] - F1[1][i]);
        u13[i] = cons_vars[2][i] - (dt/dx)*(F1[2][i+1] - F1[2][i]);
    }

    u1.push_back(u11);
    u1.push_back(u12);
    u1.push_back(u13);

    vvd F2 = CUS::get_cu_numerical_flux(u1, dt, false, dx);

    vvd u2;
    vd u21(len), u22(len), u23(len);

    for(ll i=0 ; i<len ; i++){
        u21[i] = ( 3*cons_vars[0][i] + u11[i] - (dt/dx)*(F2[0][i+1] - F2[0][i]) ) / 4;
        u22[i] = ( 3*cons_vars[1][i] + u12[i] - (dt/dx)*(F2[1][i+1] - F2[1][i]) ) / 4;
        u23[i] = ( 3*cons_vars[2][i] + u13[i] - (dt/dx)*(F2[2][i+1] - F2[2][i]) ) / 4;
    }

    u2.push_back(u21);
    u2.push_back(u22);
    u2.push_back(u23);

    vvd F3 = CUS::get_cu_numerical_flux(u2, dt, false, dx);

    for(ll i=0 ; i<len ; i++){
        cons_vars[0][i] = ( cons_vars[0][i] + 2*u21[i] - 2*(dt/dx)*(F3[0][i+1] - F3[0][i]) ) / 3;
        cons_vars[1][i] = ( cons_vars[1][i] + 2*u22[i] - 2*(dt/dx)*(F3[1][i+1] - F3[1][i]) ) / 3;
        cons_vars[2][i] = ( cons_vars[2][i] + 2*u23[i] - 2*(dt/dx)*(F3[2][i+1] - F3[2][i]) ) / 3;
    }

}




