/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/
#pragma once

#include<iostream>
#include<vector>
#include<utility>
#include<unordered_map>
#include<string>
#include<cmath>

#include "Constants.h"

#define vd vector<double>
#define vvd vector< vector<double> >

using namespace std;

namespace CTS = Constants;

namespace PrimitiveVariables{
    // Function to find the Primitive variables from Conserved variables
    vvd get_u(string problem, vvd& rho, vvd& mx);
    vvd get_v(string problem, vvd& rho, vvd& my);
    vvd get_p(string problem, vvd& rho, vvd& u, vvd& v, vvd& E);
}

/*
    Function implementations
*/
vvd PrimitiveVariables::get_u(string problem, vvd& rho, vvd& mx){
    int row = rho.size();
    int col = rho[0].size();
    
    vvd u(row, vd(col));

    for(int i=0 ; i<row ; i++){
        for(int j=0 ; j<col ; j++){
            u[i][j] = mx[i][j] / rho[i][j];
        }
    }

    return u;
}

vvd PrimitiveVariables::get_v(string problem, vvd& rho, vvd& my){
    int row = rho.size();
    int col = rho[0].size();

    vvd v(row, vd(col));

    for(int i=0 ; i<row ; i++){
        for(int j=0 ; j<col ; j++){
            v[i][j] = my[i][j] / rho[i][j];
        }
    }

    return v;
}

vvd PrimitiveVariables::get_p(string problem, vvd& rho, vvd& u, vvd& v, vvd& E){
    int row = rho.size();
    int col = rho[0].size();

    vvd p(row, vd(col));

    for(int i=0 ; i<row ; i++){
        for(int j=0 ; j<col ; j++){
            p[i][j] = ( CTS::GAMMA-1 )*( E[i][j] - 0.5*rho[i][j]*( pow(u[i][j],2) + pow(v[i][j],2) ) ) ;
        }
    }

    return p;
}

