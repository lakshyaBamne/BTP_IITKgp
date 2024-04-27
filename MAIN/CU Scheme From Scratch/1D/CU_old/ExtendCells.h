#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<utility>

#include "Constants.h"
#include "Parameters.h"

using namespace std;

namespace CTS = Constants;
namespace PRM = Parameters;

#define ll long long int
#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define pdd pair<double,double>
#define pss pair<string,string>
#define pii pair<ll,ll>

namespace ExtendCells{
    // function to extend the conserved variables
    void extend_conserved_variables(vvd &cons_vars);

    // function to extend a single vector
    void extend_vector(vd& var, string var_name);
}

void ExtendCells::extend_vector(vd &var, string var_name){
    switch(PRM::BC){
        case 1: // FREE
            var.push_back(var.back());
            var.insert(var.begin(), var[0]);

            break;

        case 2: // SOLIDWALL

            if( var_name == "MOMENTUM" ){
                var.push_back( (-1)*var.back() );
                var.insert( var.begin() , (-1)*var[0] );
            }
            else{
                var.push_back(var.back());
                var.insert(var.begin(), var[0]);
            }

            break;
        
        default:
            cout << "---ERROR--- Please enter a correct variable to extend ---" << endl;
            break;
    }
}

void ExtendCells::extend_conserved_variables(vvd &cons_vars){
    switch (PRM::BC){
        case 1: //FREE
            
            // density
            cons_vars[0].push_back(cons_vars[0].back());
            cons_vars[0].insert(cons_vars[0].begin(), cons_vars[0][0]);
            
            // momentum
            cons_vars[1].push_back(cons_vars[1].back());
            cons_vars[1].insert(cons_vars[1].begin(), cons_vars[1][0]);
            
            // energy
            cons_vars[2].push_back(cons_vars[2].back());
            cons_vars[2].insert(cons_vars[2].begin(), cons_vars[2][0]);

            break;

        case 2: // SOLIDWALL

            // density
            cons_vars[0].push_back(cons_vars[0].back());
            cons_vars[0].insert(cons_vars[0].begin(), cons_vars[0][0]);
            
            // momentum
            cons_vars[1].push_back((-1)*cons_vars[1].back());
            cons_vars[1].insert(cons_vars[1].begin(), (-1)*cons_vars[1][0]);
            
            // energy
            cons_vars[2].push_back(cons_vars[2].back());
            cons_vars[2].insert(cons_vars[2].begin(), cons_vars[2][0]);
            
            break;
        
        default:
            cout << "---ERROR--- Please enter correct boundary conditons ---" << endl;
            break;
    }
}

