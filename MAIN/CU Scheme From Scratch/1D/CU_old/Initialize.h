#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<cmath>

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
#define pii pair<ll,ll>
#define pss pair<string,string>

namespace Initialize{
    // Function to initialize the domain variables for the simulation
    void InitializeVariables(pdd &domain, pdd &time, ll &Nx);
    
    // Function to initialize the grid using the domain variables
    vd MakeGrid(pdd domain, ll Nx);

    // Function to initialize the Riemann problem variables
    vvd InitializeRiemannProblem(vd grid);
}

void Initialize::InitializeVariables(pdd &domain, pdd &time, ll &Nx){
    
    switch (PRM::PROBLEM){
        
        case 1: // MCW
            
            domain.first = 0;
            domain.second = 1;

            time.first = 0;
            time.second = 2;
            
            if( PRM::MODE == 1 ){
                Nx = 100;
            }
            else if( PRM::MODE == 2 ){
                Nx = 1000;
            }
            else{
                cout << "---ERROR--- Please enter correct mode to run ---" << endl;
            }

            break;
        
        case 2: // SCW

            domain.first = 0;
            domain.second = 1;

            time.first = 0;
            time.second = 0.012;
            
            if( PRM::MODE == 1 ){
                Nx = 200;
            }
            else if( PRM::MODE == 2 ){
                Nx = 1000;
            }
            else{
                cout << "---ERROR--- Please enter correct mode to run ---" << endl;
            }

            break;
        
        case 3: // LAX

            domain.first = -5;
            domain.second = 5;

            time.first = 0;
            time.second = 1.3;
            
            if( PRM::MODE == 1 ){
                Nx = 200;
            }
            else if( PRM::MODE == 2 ){
                Nx = 1000;
            }
            else{
                cout << "---ERROR--- Please enter correct mode to run ---" << endl;
            }

            break;
        
        case 4: // SOS

            break;
        
        case 5: // BLW

            domain.first = 0;
            domain.second = 1;

            time.first = 0;
            time.second = 0.038;
            
            if( PRM::MODE == 1 ){
                Nx = 1000;
            }
            else if( PRM::MODE == 2 ){
                Nx = 4000;
            }
            else{
                cout << "---ERROR--- Please enter correct mode to run ---" << endl;
            }

            break;

        default:
            cout << "---ERROR--- Please enter the correct problem ---" << endl;
            break;

    }

}

vd Initialize::MakeGrid(pdd domain, ll Nx){
    vd grid;

    grid.push_back(domain.first);

    for(ll i=0 ; i<Nx ; i++){
        grid.push_back( grid.back() + ( (domain.second-domain.first)/Nx ) );
    }

    return grid;
}

vvd Initialize::InitializeRiemannProblem(vd grid){
    vvd cons_vars; // 2D vector to store the conserved variables
    vd rho, m, E; // 1D vectors for the individual conserved variables
    double u, p;

    switch(PRM::PROBLEM){
        case 1: // MCW
            u = 0.1;
            p = 1;

            for(auto i : grid){
                if( i < 0.3 ){
                    rho.push_back( 1.4 );
                }
                else{
                    rho.push_back( 1.0 );
                }

                m.push_back( rho.back()*u );
                E.push_back( ( p/(CTS::GAMMA-1) ) + (rho.back()*(pow(u,2)/2)) );
            }

            break;
        
        case 2: // SCW
            u = -19.59745;

            for(auto i : grid){
                if( i < 0.8 ){
                    rho.push_back( 1.0 );
                    p = 1000.0;
                }
                else{
                    rho.push_back( 1.0 );
                    p = 0.01;
                }

                m.push_back( rho.back()*u );
                E.push_back( ( p/(CTS::GAMMA-1) ) + (rho.back()*(pow(u,2)/2)) );
            }
            
            break;

        case 3: // LAX
            for(auto i : grid){
                if( i < 0.5 ){
                    rho.push_back( 0.445 );
                    u = 0.698;
                    p = 3.528;
                }
                else{
                    rho.push_back( 0.5 );
                    u = 0.0;
                    p = 0.571;
                }

                m.push_back( rho.back()*u );
                E.push_back( ( p/(CTS::GAMMA-1) ) + (rho.back()*(pow(u,2)/2)) );
            }

            break;

        case 4: // SOS
            
            break;

        case 5: // BLW
            
            for(auto i : grid){
                if( i < 0.1 ){
                    rho.push_back( 1.0 );
                    u = 0.0;
                    p = 1000.0;
                }
                else if( i>=0.1 && i<=0.9 ){
                    rho.push_back( 1.0 );
                    u = 0.0;
                    p = 0.01;
                }
                else{
                    rho.push_back( 1.0 );
                    u = 0.0;
                    p = 100.0;
                }

                m.push_back( rho.back()*u );
                E.push_back( ( p/(CTS::GAMMA-1) ) + (rho.back()*(pow(u,2)/2)) );
            }

            break;
        
        default:
            cout << "---ERROR--- Please enter correct problem ---" << endl;
            break;
    }

    cons_vars.push_back(rho);
    cons_vars.push_back(m);
    cons_vars.push_back(E);

    return cons_vars;
}


