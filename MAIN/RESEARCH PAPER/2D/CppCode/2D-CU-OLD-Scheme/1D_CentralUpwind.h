#pragma once

#include<iostream>
#include<vector>
#include<utility>
#include<string>
#include<cmath>
#include<sstream>
#include<fstream>

#include "Constants.h"
#include "Utility.h"
#include "InitializeRiemannProblem.h"
#include "ExtendCells.h"
#include "CUNumericalFlux.h"
#include "SSPRK.h"
#include "OutputResult.h"

namespace IRP = InitRiemannProblem;
namespace UTL = Utility;
namespace EXC = ExtendCells;
namespace CTS = Constants;
namespace CUF = CUNumericalFlux;
namespace SRK = SSPRK;
namespace OPF = OutputResult;

#define ll long long int

using namespace std;

class GetInput{

public: // variables
    pair<double,double> domainX, domainY;
    pair<double,double> time;
    ll Nx, Ny;
    pair<string,string> initial_conditions;

    double dt;
    double t;
    double dx, dy;

public: // methods
    GetInput(){
        // show to the user information about the lubrary and it's basic usage
        show_info();

        // take input for the starting state
        cout << "|---------------------------------------INPUT----------------------------------------|" << endl;
        cout << "+------------------------------------------------------------------------------------+" << endl;

        cout << "DomainX start : ";
        cin >> domainX.first;

        cout << "DomainX end : ";
        cin >> domainX.second;

        cout << "DomainY start : ";
        cin >> domainY.first;

        cout << "DomainY end : ";
        cin >> domainY.second;

        cout << "Grid points in domainX : ";
        cin >> Nx;

        cout << "Grid points in domainY : ";
        cin >> Ny;

        cout << "Initial time : ";
        cin >> time.first;

        cout << "Final time : ";
        cin >> time.second;

        cout << "Enter the Riemann Problem : ";
        cin >> initial_conditions.first;

        cout << "Enter the Boundary Conditions : ";
        cin >> initial_conditions.second;

        cout << "|------------------------------------------------------------------------------------|" << endl;

        t = time.first;

        dx = ( domainX.second - domainX.first ) / Nx;
        dy = ( domainY.second - domainY.first ) / Ny;

        stringstream ss1, ss2;
        ss1 << Nx << " " << Ny;

        UTL::export_string(initial_conditions.first);
        UTL::export_string(ss1.str());

        ss2 << domainX.first << " " << domainX.second << " " << domainY.first << " " << domainY.second;
        UTL::export_string(ss2.str());
    }

    // function to show information about the problems solvable by the system
    void show_info(){
        cout << "+------------------------------------------------------------------------------------+" << endl;
        cout << "|                       2-Dimensional Central Upwind Scheme                          |" << endl;
        cout << "+------------------------------------------------------------------------------------+" << endl;
        cout << "|                  Riemann problems solvable using this library                      |" << endl;
        cout << "+------------------------------------------------------------------------------------+" << endl;

        cout << "-> [MCW] Slowly moving isolated contact discontinuity" << endl;
        cout << "-> [2DR-CFG3] 2 Dimensional Riemann Problem (Configuration-3)" << endl;
        cout << "-> [EXP] Explosion Problem" << endl;

        cout << "+------------------------------------------------------------------------------------+" << endl;
        cout << "|                            Boundary Conditions supported                           |" << endl;
        cout << "+------------------------------------------------------------------------------------+" << endl;
        
        cout << "-> [FREE] Neumann Boundary Conditions" << endl;
        cout << "-> [RFREE] Reflect Free Boundary Conditions" << endl;

        cout << "+------------------------------------------------------------------------------------+" << endl;
    }

    //! Method to run the old CU Scheme for 2D using the same code for 1D
    void run_cu_old_partial(string result){
        vector< vector< vector<double> > > cons_vars = IRP::get_conserved_variables(initial_conditions, domainX, domainY, Nx, Ny, result);

        double amax=0;
        double bmax=0;    

        string file_name = result + "/density0.txt";
        OPF::write_density(cons_vars[0], file_name);

        while( t < time.second ){

            // Now we want to run the CU Scheme for this 2D grid using the 1D code
            vector< vector< vector<double> > > F = get_flux_x(cons_vars, amax);
            vector< vector< vector<double> > > G = get_flux_y(cons_vars, bmax);

            // update time step
            dt = CTS::CFL*min(dx/amax, dy/bmax);

            // update the conserved variables
            for(int u=0 ; u<4 ; u++){
                for(int i=0 ; i<Nx ; i++){
                    for(int j=0 ; j<Ny ; j++){
                        cons_vars[u][i][j] = cons_vars[u][i][j] - (dt/dx)*( F[j][u][i+1] - F[j][u][i] ) - (dt/dy)*( G[i][u][j+1] - G[i][u][j] );
                    }
                }
            }

            t += dt;
            cout << "t=" << t << " | dt=" << dt << endl;
        }

        file_name = result + "/density1.txt";
        OPF::write_density(cons_vars[0], file_name);

    }
    
    //! Method to calculate the CU Numerical Flux in x-direction
    vector< vector< vector<double> > > get_flux_x(vector< vector< vector<double> > >& U, double& amax){
        vector< vector< vector<double> > > F;

        double amax_global;
        double amax_local;


        // we want the CU Flux for all the rows in the domain
        for(int j=0 ; j<Ny ; j++){
            vector< vector<double> > cons_vars(4, vector<double>(Nx, 0));

            for(int u=0 ; u<4 ; u++){
                for(int i=0 ; i<Nx ; i++){
                    cons_vars[u][i] = U[u][i][j];
                }
            }

            EXC::extend_cells_horizontal(initial_conditions.second, cons_vars);

            F.push_back( CUF::get_cu_flux_horizontal(cons_vars, initial_conditions, dt, dx, t, time, amax_local) );

            amax_global = max(amax_global, amax_local);
        }

        // update the amax value
        amax = amax_global;

        return F;

    }

    //! Method to calculate the CU Numerical Flux in y-direction
    vector< vector< vector<double> > > get_flux_y(vector< vector< vector<double> > >& U, double& bmax){

        vector< vector< vector<double> > > G;

        double bmax_global;
        double bmax_local;

        // we want the CU Flux for all the rows in the domain
        for(int i=0 ; i<Nx ; i++){
            vector< vector<double> > cons_vars(4, vector<double>(Ny, 0));

            for(int u=0 ; u<4 ; u++){
                for(int j=0 ; j<Ny ; j++){
                    cons_vars[u][j] = U[u][i][j];
                }
            }

            EXC::extend_cells_vertical(initial_conditions.second, cons_vars);

            G.push_back( CUF::get_cu_flux_vertical(cons_vars, initial_conditions, dt, dx, t, time, bmax_local) );

            bmax_global = max(bmax_global, bmax_local);
        }

        // update the amax value
        bmax = bmax_global;

        return G;

    }
    
};

