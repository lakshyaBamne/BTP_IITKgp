/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>
#include<sstream>
#include<fstream>

#include "ExtendCells.h"
#include "OutputResult.h"
#include "Constants.h"
#include "PrimitiveVariables.h"
#include "Utility.h"
#include "Flux.h"

#define ll long long int
#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define vvvvd vector< vector< vector< vector<double> > > >
#define umap_ss unordered_map<string,string>

using namespace std;

namespace EXC = ExtendCells;
namespace OPR = OutputResult;
namespace CTS = Constants;
namespace PRV = PrimitiveVariables;
namespace UTL = Utility;
namespace FLX = Flux;

// class RiemannProblem(RP) to initialize the 2D riemann problem for a given problem
class RP{
public: // everything is public

    // required later for various purposes
    string mode;
    string problem;

    pair<double,double> domX;
    pair<double,double> domY;

    pair<double,double> time;

    int Nx;
    int Ny;

    double dx;
    double dy;

    double t;
    double dt; // time discretization step updated using CFL conditions

    // grid vectors are initialized after initializing the variables
    vector< double > gridX;
    vector< double > gridY;

    // Variable matrices to store information
    vector< vector<double> > rho;   
    vector< vector<double> > mx;
    vector< vector<double> > my;
    vector< vector<double> > E;

    string bc;

    /*
        Constructor takes the problem string as input to initialize the variables
    
        -> [MCW] Moving Contact Wave
        -> [2DR] 2-Dimensional Riemann Problem
        -> [EXP] Explosion Problem
        -> [IMP1, IMP2] Implosion Problem
        -> [KHI1, KHI2, KHI3] Kelvin Helmholtz Instability
        -> [RTI] Raleigh Taylor Instability
    */
    RP(string mode, string problem) : mode(mode), problem(problem){

        cout << "---LOG--- Constructor running ---" << endl;

        // initialize the environment and set up initial variables and matrices
        initialize_environment();        

        cout << "---LOG--- environment write successful ---" << endl;
    }

    // Function to initialize environment variables and set up initial matrices
    void initialize_environment(){
        if( mode == "CU" ){
            // CU plot is the plot for the Numerical Scheme itself
            if( problem == "MCW" ){
                domX.first = -0.2;
                domX.second = 0.2;

                domY.first = 0;
                domY.second = 0.8;

                dx = (double)1/250;
                dy = (double)1/250;

                time.first = 0;
                time.second = 2;

                t = time.first;

                bc = "FREE";
            }
            else if( problem == "2DR" ){
                domX.first = 0;
                domX.second = 1.2;

                domY.first = 0;
                domY.second = 1.2;

                dx = (double)3/500;
                dy = (double)3/500;

                time.first = 0;
                time.second = 1;

                t = time.first;

                bc = "FREE";
            }
            else if( problem == "EXP" ){
                domX.first = 0;
                domX.second = 1.5;

                domY.first = 0;
                domY.second = 1.5;

                dx = (double)3/800;
                dy = (double)3/800;

                time.first = 0;
                time.second = 3.2;

                t = time.first;

                bc = "RFREE";
            }
            else if( problem == "IMP" ){
                domX.first = 0;
                domX.second = 0.3;

                domY.first = 0;
                domY.second = 0.3;

                dx = (double)3/6000;
                dy = (double)3/6000;

                // dx = (double)3/600;
                // dy = (double)3/800;

                time.first = 0;
                time.second = 2.5;

                t = time.first;

                bc = "SOLIDWALL";
            }
            else if( problem == "KHI" ){
                domX.first = -0.5;
                domX.second = 0.5;

                domY.first = -0.5;
                domY.second = 0.5;

                dx = (double)3/1024;
                dy = (double)3/1024;

                time.first = 0;
                time.second = 1; // also evaluate results at t=2.5 and t=4

                t = time.first;

                bc = "PERIODIC";
            }
            else if( problem == "RTI" ){
                domX.first = 0;
                domX.second = 0.25;

                domY.first = 0;
                domY.second = 1.0;

                dx = (double)1/1024;
                dy = (double)1/1024;

                time.first = 0;
                time.second = 2.95; // also evaluate results at t=2.5 and t=4

                t = time.first;

                bc = "DIRICHLET";
            }
            else{
                cout << "---ERROR--- Please enter correct problem ---" << endl;
            }
        }
        else if( mode == "REF" ){
            // reference is taken as 4 times finer grid
            if( problem == "MCW" ){
                domX.first = -0.2;
                domX.second = 0.2;

                domY.first = 0;
                domY.second = 0.8;

                dx = (double)1/1000;
                dy = (double)1/1000;

                time.first = 0;
                time.second = 2;

                t = time.first;

                bc = "FREE";
            }
            else if( problem == "2DR" ){
                domX.first = 0;
                domX.second = 1.2;

                domY.first = 0;
                domY.second = 1.2;

                dx = (double)3/4000;
                dy = (double)3/4000;

                time.first = 0;
                time.second = 1;

                t = time.first;

                bc = "FREE";
            }
            else if( problem == "EXP" ){
                domX.first = 0;
                domX.second = 1.5;

                domY.first = 0;
                domY.second = 1.5;

                dx = (double)3/3200;
                dy = (double)3/3200;

                time.first = 0;
                time.second = 3.2;

                t = time.first;

                bc = "RFREE";
            }
            else if( problem == "IMP" ){
                domX.first = 0;
                domX.second = 0.3;

                domY.first = 0;
                domY.second = 0.3;

                dx = (double)3/24000;
                dy = (double)3/24000;

                time.first = 0;
                time.second = 2.5;

                t = time.first;

                bc = "SOLIDWALL";
            }
            else if( problem == "KHI" ){
                domX.first = -0.5;
                domX.second = 0.5;

                domY.first = -0.5;
                domY.second = 0.5;

                dx = (double)3/4096;
                dy = (double)3/4096;

                time.first = 0;
                time.second = 1; // also evaluate results at t=2.5 and t=4

                t = time.first;

                bc = "PERIODIC";
            }
            else if( problem == "RTI" ){
                domX.first = 0;
                domX.second = 0.25;

                domY.first = 0;
                domY.second = 1.0;

                dx = (double)1/4096;
                dy = (double)1/4096;

                time.first = 0;
                time.second = 2.95; // also evaluate results at t=2.5 and t=4

                t = time.first;

                bc = "DIRICHLET";
            }
            else{
                cout << "---ERROR--- Please enter correct problem ---" << endl;
            }
        }
        else{
            cout << "---ERROR--- Please enter correct mode for initialization ---" << endl;
        }

        // calculating the number of points in the grid using the cell volumes
        Nx = (int)( (domX.second-domX.first) / dx );
        Ny = (int)( (domY.second-domY.first) / dy );

        cout << "---LOG--- Initialized Riemann problem variables ---" << endl;

        // initialize grid vectors
        make_grids();

        cout << "---LOG--- Initialized computational grid ---" << endl;

        OPR::write_grids(mode, gridX, gridY);

        // initialize variable matrices
        init_vars();

        cout << "---LOG--- Initialized Variable matrices ---" << endl;

        // export environment
        export_env(mode);
    }

    // function to initialize the computational grid
    void make_grids(){
        // gridX
        for(int i=0 ; i<Nx ; i++){
            gridX.push_back( domX.first + dx*( i + 0.5 ) );
        }

        // gridY
        for(int i=0 ; i<Ny ; i++){
            gridY.push_back( domY.first + dy*( i + 0.5 ) );
        }
    }

    // function to initialize the variable grids (rho, mx, my, E)
    void init_vars(){

        if( problem == "MCW" ){
            
            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    bool c1 = i>-0.1 && i<0.1 && j>0 && j<0.02;
                    bool c2 = i>-0.02 && i<0.02 && j>0.02 && j<0.1;
                    bool c3 = pow(i+0.02,2) + pow(j-0.02,2) < pow(0.08,2);
                    bool c4 = pow(i-0.02,2) + pow(j-0.02,2) < pow(0.08,2);

                    if( c1 || c2 || c3 || c4 ){
                        rho1.push_back( 1.4 );
                    }
                    else{
                        rho1.push_back( 1.0 );
                    }

                    u = 0.0;
                    v = 0.2;
                    p = 1;

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }
        }
        else if( problem == "2DR" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( i>1 && j>1 ){
                        rho1.push_back( 1.5 );
                        u = 0.0;
                        v = 0.0;
                        p = 1.5;
                    }
                    else if( i<1 && j>1 ){
                        rho1.push_back( 0.5323 );
                        u = 1.206;
                        v = 0.0;
                        p = 0.3;
                    }
                    else if( i<1 && j<1 ){
                        rho1.push_back( 0.138 );
                        u = 1.206;
                        v = 1.206;
                        p = 0.029;
                    }
                    else if( i>1 && j<1 ){
                        rho1.push_back( 0.5323 );
                        u = 0.0;
                        v = 1.206;
                        p = 0.3;
                    }

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }
        else if( problem == "EXP" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( pow(i,2)+pow(j,2) < 0.16 ){
                        rho1.push_back(1);
                        u = 0;
                        v = 0;
                        p = 1;
                    }
                    else{
                        rho1.push_back(0.125);
                        u = 0;
                        v = 0;
                        p = 0.1;
                    }

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }
        else if( problem == "IMP" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( abs(i) + abs(j) < 0.15 ){
                        rho1.push_back(0.125);
                        u = 0;
                        v = 0;
                        p = 0.14;
                    }
                    else{
                        rho1.push_back(1);
                        u = 0;
                        v = 0;
                        p = 1;
                    }

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }
        else if( problem == "KHI" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( j>=-0.5 && j<-0.25 ){
                        rho1.push_back(1);
                        u = -0.5 + 0.5*expf( (j+0.25)/CTS::SMOOTHP );
                    }
                    else if( j>=-0.25 && j<0 ){
                        rho1.push_back(2);
                        u = 0.5 - 0.5*expf( (-j-0.25)/CTS::SMOOTHP );
                    }
                    else if( j>=0 && j<0.25 ){
                        rho1.push_back(2);
                        u = 0.5 - 0.5*expf( (j-0.25)/CTS::SMOOTHP );
                    }
                    else if( j>=0.25 && j<=0.5 ){
                        rho1.push_back(1);
                        u = -0.5 + 0.5*expf( (0.25-j)/CTS::SMOOTHP );
                    }
                    
                    v = 0.01*sin(4*CTS::PI*i);
                    p = 1.5;

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }
        else if( problem == "RTI" ){
            
        }   
        else{
            cout << "---ERROR--- Please enter correct problem ---" << endl;
        }

        cout << "---LOG--- Initialized the variables ---" << endl;

        // extend cells
        EXC::extend_matrix(rho, bc, "Density");
        EXC::extend_matrix(mx, bc, "MomentumX");
        EXC::extend_matrix(my, bc, "MomentumY");
        EXC::extend_matrix(E, bc, "Energy");
        
        cout << "---LOG--- Extended cells successfully ---" << endl;

    }

    // function to export the environment variables in a file -> Environment.txt
    void export_env(string mode){
        string file_name;

        if( mode == "CU" ){
            file_name = "env1/Environment.txt";
        }
        else{
            file_name = "env2/Environment.txt";
        }
        
        stringstream ss;
        
        ofstream fout;
        fout.open(file_name , ios::app);

        ss << mode << endl;
        ss << problem << endl;
        ss << domX.first << " " << domX.second << endl;
        ss << domY.first << " " << domY.second << endl;
        ss << time.first << " " << time.second << endl;
        ss << Nx << endl;
        ss << Ny << endl;

        fout << ss.str() << endl;

        fout.close();
    }

    // function to run the CU Scheme
    // and store data only for initial and final iterations
    void RunCU_partial(){

        //! write initial values for plotting
        OPR::write_matrix(mode, "Density", rho);

        while( t < time.second ){
            cout << "t = " << t << " | dt = " << dt << endl;
            
            // function to update the conservative variables internally
            // also update dt for next iteration
            get_next_vars();

            // update new time step calculated using CFL conditions
            t += dt;

        }

        //! write final values for plotting
        OPR::write_matrix(mode, "Density", rho);
    }

    /*
        Functions to Run the CU Scheme
    */
    void get_next_vars(){
        // using class variables as input
        vvvvd cu_flux = get_cu_flux(rho, mx, my, E, true);
        
        vvvd F = cu_flux[0];
        vvvd G = cu_flux[1];

        double LAMBDA = dt/dx;
        double MU = dt/dy;

        int row = rho.size();
        int col = rho[0].size();

        // STAGE-1
        vvd rho1(row-2, vd(col-2));
        vvd mx1(row-2, vd(col-2));
        vvd my1(row-2, vd(col-2));
        vvd E1(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                rho1[i-1][j-1] = rho[i][j] - ( LAMBDA*( F[0][i][j] - F[0][i][j-1] ) + MU*( G[0][i][j] - G[0][i-1][j] ) );
                mx1[i-1][j-1] = mx[i][j] - ( LAMBDA*( F[1][i][j] - F[1][i][j-1] ) + MU*( G[1][i][j] - G[1][i-1][j] ) );
                my1[i-1][j-1] = my[i][j] - ( LAMBDA*( F[2][i][j] - F[2][i][j-1] ) + MU*( G[2][i][j] - G[2][i-1][j] ) );
                E1[i-1][j-1] = E[i][j] - ( LAMBDA*( F[3][i][j] - F[3][i][j-1] ) + MU*( G[3][i][j] - G[3][i-1][j] ) );
            }
        }

        EXC::extend_matrix(rho1, bc, "Density");
        EXC::extend_matrix(mx1, bc, "MomentumX");
        EXC::extend_matrix(my1, bc, "MomentumY");
        EXC::extend_matrix(E1, bc, "Energy");

        // STAGE-2
        vvvvd cu_flux_1 = get_cu_flux(rho1, mx1, my1, E1, false);

        vvvd F1 = cu_flux_1[0];
        vvvd G1 = cu_flux_1[1];

        vvd rho2(row-2, vd(col-2));
        vvd mx2(row-2, vd(col-2));
        vvd my2(row-2, vd(col-2));
        vvd E2(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                rho2[i-1][j-1] = ( 3*rho[i][j] + rho1[i][j] - ( LAMBDA*( F1[0][i][j] - F1[0][i][j-1] ) + MU*( G1[0][i][j] - G1[0][i-1][j] ) ) ) / 4;
                mx2[i-1][j-1] = ( 3*mx[i][j] + mx1[i][j] - ( LAMBDA*( F1[1][i][j] - F1[1][i][j-1] ) + MU*( G1[1][i][j] - G1[1][i-1][j] ) ) ) / 4;
                my2[i-1][j-1] = ( 3*my[i][j] + my1[i][j] - ( LAMBDA*( F1[2][i][j] - F1[2][i][j-1] ) + MU*( G1[2][i][j] - G1[2][i-1][j] ) ) ) / 4;
                E2[i-1][j-1] = ( 3*E[i][j] + E1[i][j] - ( LAMBDA*( F1[3][i][j] - F1[3][i][j-1] ) + MU*( G1[3][i][j] - G1[3][i-1][j] ) ) ) / 4;
            }
        }

        EXC::extend_matrix(rho2, bc, "Density");
        EXC::extend_matrix(mx2, bc, "MomentumX");
        EXC::extend_matrix(my2, bc, "MomentumY");
        EXC::extend_matrix(E2, bc, "Energy");

        // STAGE-3
        vvvvd cu_flux_new = get_cu_flux(rho2, mx2, my2, E2, false);

        vvvd FN = cu_flux_new[0];
        vvvd GN = cu_flux_new[1];

        vvd rhon(row-2, vd(col-2));
        vvd mxn(row-2, vd(col-2));
        vvd myn(row-2, vd(col-2));
        vvd En(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                rhon[i-1][j-1] = ( rho[i][j] + 2*rho2[i][j] - 2*( LAMBDA*( FN[0][i][j] - FN[0][i][j-1] ) + MU*( GN[0][i][j] - GN[0][i-1][j] ) ) ) / 3;
                mxn[i-1][j-1] = ( mx[i][j] + 2*mx2[i][j] - 2*( LAMBDA*( FN[1][i][j] - FN[1][i][j-1] ) + MU*( GN[1][i][j] - GN[1][i-1][j] ) ) ) / 3;
                myn[i-1][j-1] = ( my[i][j] + 2*my2[i][j] - 2*( LAMBDA*( FN[2][i][j] - FN[2][i][j-1] ) + MU*( GN[2][i][j] - GN[2][i-1][j] ) ) ) / 3;
                En[i-1][j-1] = ( E[i][j] + 2*E2[i][j] - 2*( LAMBDA*( FN[3][i][j] - FN[3][i][j-1] ) + MU*( GN[3][i][j] - GN[3][i-1][j] ) ) ) / 3;
            }
        }

        EXC::extend_matrix(rhon, bc, "Density");
        EXC::extend_matrix(mxn, bc, "MomentumX");
        EXC::extend_matrix(myn, bc, "MomentumY");
        EXC::extend_matrix(En, bc, "Energy");

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col ; j++){
                rho[i][j] = rhon[i][j];
                mx[i][j] = mxn[i][j];
                my[i][j] = myn[i][j];
                E[i][j] = En[i][j];
            }
        }

    }

    // function to get the CU flux final for any given 2D grid of Conserved variables
    // no arguments required because the grids are already class variables accessible to all functions
    //! returns the flux matrices F and G which are then used in SSPRK scheme to construct next iteration
    vvvvd get_cu_flux(vvd& rho, vvd& mx, vvd& my, vvd& E, bool update_time){
        int row = rho.size();
        int col = rho[0].size();
        
        //! Piecewise Linear Reconstructions
        vvvd rho_plr = get_plr(rho, "Density");
        vvvd mx_plr = get_plr(mx, "MomentumX");
        vvvd my_plr = get_plr(my, "MomentumY");
        vvvd E_plr = get_plr(E, "Energy");

        //! Local Speeds of Propagations at the Interfaces
        vvvd lsp_x = get_lsp_a(rho_plr, mx_plr, my_plr, E_plr);
        vvvd lsp_y = get_lsp_b(rho_plr, mx_plr, my_plr, E_plr);

        //! Flux Vector values
        // vvvd fflux_E = FLX::fflux(problem, rho_plr[2], mx_plr[2], my_plr[2], E_plr[2]);
        // vvvd fflux_W = FLX::fflux(problem, rho_plr[3], mx_plr[3], my_plr[3], E_plr[3]);        
        // vvvd gflux_N = FLX::gflux(problem, rho_plr[0], mx_plr[0], my_plr[0], E_plr[0]);
        // vvvd gflux_S = FLX::gflux(problem, rho_plr[1], mx_plr[1], my_plr[1], E_plr[1]);

        //! Update the time step
        if( update_time ){
            new_dt(lsp_x, lsp_y);
        }

        //! Get CU Numerical Flux in Horizontal direction        
        vvvd F;

        vvd F1(row-2, vd(col-1));
        vvd F2(row-2, vd(col-1));
        vvd F3(row-2, vd(col-1));
        vvd F4(row-2, vd(col-1));

        double uE, uW, vE, vW, pE, pW;
        double f1E, f2E, f3E, f4E;
        double f1W, f2W, f3W, f4W;

        for(int i=1 ; i<row-2 ; i++){
            for(int j=0 ; j<col-1 ; j++){
                uE = mx_plr[2][i][j] / rho_plr[2][i][j];
                uW = mx_plr[3][i][j+1]/rho_plr[3][i][j+1];
                vE = my_plr[2][i][j] / rho_plr[2][i][j];
                vW = my_plr[3][i][j+1] / rho_plr[3][i][j+1];
                pE = (CTS::GAMMA-1)*( E_plr[2][i][j] - 0.5*rho_plr[2][i][j]*( pow(uE,2) + pow(vE,2) ) );
                pW = (CTS::GAMMA-1)*( E_plr[3][i][j+1] - 0.5*rho_plr[3][i][j+1]*( pow(uW,2) + pow(vW,2) ) );

                f1E = mx_plr[2][i][j];
                f2E = mx_plr[2][i][j]*uE + pE;
                f3E = mx_plr[2][i][j]*vE;
                f4E = uE*( E_plr[2][i][j] + pE );

                f1W = mx_plr[3][i][j+1];
                f2W = mx_plr[3][i][j+1]*uW + pW;
                f3W = mx_plr[3][i][j+1]*vW;
                f4W = uW*( E_plr[3][i][j+1] + pW );

                double diff = lsp_x[0][i][j] - lsp_x[1][i][j];
                double prod = lsp_x[0][i][j]*lsp_x[1][i][j];
                
                if( diff > CTS::EPSILON ){

                    // F1[i][j] = ( lsp_x[0][i][j]*fflux_E[0][i][j] - lsp_x[1][i][j]*fflux_W[0][i][j+1] + prod*( rho_plr[3][i][j+1] - rho_plr[2][i][j] ) ) / (diff);
                    // F2[i][j] = ( lsp_x[0][i][j]*fflux_E[1][i][j] - lsp_x[1][i][j]*fflux_W[1][i][j+1] + prod*( mx_plr[3][i][j+1] - mx_plr[2][i][j] ) ) / (diff);
                    // F3[i][j] = ( lsp_x[0][i][j]*fflux_E[2][i][j] - lsp_x[1][i][j]*fflux_W[2][i][j+1] + prod*( my_plr[3][i][j+1] - my_plr[2][i][j] ) ) / (diff);
                    // F4[i][j] = ( lsp_x[0][i][j]*fflux_E[3][i][j] - lsp_x[1][i][j]*fflux_W[3][i][j+1] + prod*( E_plr[3][i][j+1] - E_plr[2][i][j] ) ) / (diff);

                    F1[i][j] = ( lsp_x[0][i][j]*f1E - lsp_x[1][i][j]*f1W + prod*( rho_plr[3][i][j+1] - rho_plr[2][i][j] ) ) / (diff);
                    F2[i][j] = ( lsp_x[0][i][j]*f2E - lsp_x[1][i][j]*f2W + prod*( mx_plr[3][i][j+1] - mx_plr[2][i][j] ) ) / (diff);
                    F3[i][j] = ( lsp_x[0][i][j]*f3E - lsp_x[1][i][j]*f3W + prod*( my_plr[3][i][j+1] - my_plr[2][i][j] ) ) / (diff);
                    F4[i][j] = ( lsp_x[0][i][j]*f4E - lsp_x[1][i][j]*f4W + prod*( E_plr[3][i][j+1] - E_plr[2][i][j] ) ) / (diff);

                }
                else{

                    // F1[i][j] = 0.5*(fflux_E[0][i][j]+fflux_W[0][i][j]);
                    // F2[i][j] = 0.5*(fflux_E[1][i][j]+fflux_W[1][i][j]);
                    // F3[i][j] = 0.5*(fflux_E[2][i][j]+fflux_W[2][i][j]);
                    // F4[i][j] = 0.5*(fflux_E[3][i][j]+fflux_W[3][i][j]);

                    F1[i][j] = 0.5*(f1E + f1W);
                    F2[i][j] = 0.5*(f2E + f2W);
                    F3[i][j] = 0.5*(f3E + f3W);
                    F4[i][j] = 0.5*(f4E + f4W);

                }

            }
        }

        F.push_back(F1);
        F.push_back(F2);
        F.push_back(F3);
        F.push_back(F4);

        //! Get CU Numerical Flux in Vertical direction
        vvvd G;

        vvd G1(row-1, vd(col-2));
        vvd G2(row-1, vd(col-2));
        vvd G3(row-1, vd(col-2));
        vvd G4(row-1, vd(col-2));

        double uN, uS, vN, vS, pN, pS;
        double g1N, g2N, g3N, g4N;
        double g1S, g2S, g3S, g4S;

        for(int i=0 ; i<row-1 ; i++){
            for(int j=1 ; j<col-2 ; j++){

                uN = mx_plr[0][i][j] / rho_plr[0][i][j];
                uS = mx_plr[1][i+1][j] / rho_plr[1][i+1][j];
                vN = my_plr[0][i][j] / rho_plr[0][i][j];
                vS = my_plr[1][i+1][j] / rho_plr[1][i+1][j];
                pN = (CTS::GAMMA-1)*( E_plr[0][i][j] - 0.5*rho_plr[0][i][j]*( pow(uN,2) + pow(vN,2) ) );
                pS = (CTS::GAMMA-1)*( E_plr[1][i+1][j] - 0.5*rho_plr[1][i+1][j]*( pow(uS,2) + pow(vS,2) ) );

                g1N = my_plr[0][i][j];
                g2N = my_plr[0][i][j]*uN;
                g3N = my_plr[0][i][j]*vN + pN;
                g4N = vN*( E_plr[0][i][j] + pN );

                g1S = my_plr[1][i+1][j];
                g2S = my_plr[1][i+1][j]*uS;
                g2S = my_plr[1][i+1][j]*vS + pS;
                g4S = vS*( E_plr[1][i+1][j] + pS );

                double diff = lsp_y[0][i][j] - lsp_y[1][i][j];
                double prod = lsp_y[0][i][j]*lsp_y[1][i][j];
                
                if( diff > CTS::EPSILON ){
                    // use ADT

                    // G1[i][j] = ( lsp_y[0][i][j]*gflux_N[0][i][j] - lsp_y[1][i][j]*gflux_S[0][i+1][j] + prod*( rho_plr[1][i+1][j] - rho_plr[0][i][j] ) ) / (diff);
                    // G2[i][j] = ( lsp_y[0][i][j]*gflux_N[1][i][j] - lsp_y[1][i][j]*gflux_S[1][i+1][j] + prod*( mx_plr[1][i+1][j] - mx_plr[0][i][j] ) ) / (diff);
                    // G3[i][j] = ( lsp_y[0][i][j]*gflux_N[2][i][j] - lsp_y[1][i][j]*gflux_S[2][i+1][j] + prod*( my_plr[1][i+1][j] - my_plr[0][i][j] ) ) / (diff);
                    // G4[i][j] = ( lsp_y[0][i][j]*gflux_N[3][i][j] - lsp_y[1][i][j]*gflux_S[3][i+1][j] + prod*( E_plr[1][i+1][j] - E_plr[0][i][j] ) ) / (diff);

                    G1[i][j] = ( lsp_y[0][i][j]*g1N - lsp_y[1][i][j]*g1S + prod*( rho_plr[1][i+1][j] - rho_plr[0][i][j] ) ) / (diff);
                    G2[i][j] = ( lsp_y[0][i][j]*g2N - lsp_y[1][i][j]*g2S + prod*( mx_plr[1][i+1][j] - mx_plr[0][i][j] ) ) / (diff);
                    G3[i][j] = ( lsp_y[0][i][j]*g3N - lsp_y[1][i][j]*g3S + prod*( my_plr[1][i+1][j] - my_plr[0][i][j] ) ) / (diff);
                    G4[i][j] = ( lsp_y[0][i][j]*g4N - lsp_y[1][i][j]*g4S + prod*( E_plr[1][i+1][j] - E_plr[0][i][j] ) ) / (diff);

                }
                else{
                    // use normal averaging

                    // G1[i][j] = 0.5*(gflux_N[0][i][j]+gflux_S[0][i][j]);
                    // G2[i][j] = 0.5*(gflux_N[1][i][j]+gflux_S[1][i][j]);
                    // G3[i][j] = 0.5*(gflux_N[2][i][j]+gflux_S[2][i][j]);
                    // G4[i][j] = 0.5*(gflux_N[3][i][j]+gflux_S[3][i][j]);

                    G1[i][j] = 0.5*(g1N + g1S);
                    G2[i][j] = 0.5*(g2N + g2S);
                    G3[i][j] = 0.5*(g3N + g3S);
                    G4[i][j] = 0.5*(g4N + g4S);

                }

            }
        }

        G.push_back(G1);
        G.push_back(G2);
        G.push_back(G3);
        G.push_back(G4);

        return {F,G};
    }

    // Function to calculate new time step using the CFL conditions
    void new_dt(vvvd& lsp_x, vvvd& lsp_y){
        double amax=0; 
        double bmax=0;

        // both lsp_x and lsp_y have the same dimensions but different relevant values
        int row = lsp_x[0].size();
        int col = lsp_x[0][0].size();

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){
                amax = max(
                    amax,
                    max(
                        lsp_x[0][i][j], -1*lsp_x[1][j][j]
                    )
                );
            }
        }

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){
                bmax = max(
                    bmax,
                    max(
                        lsp_y[0][i][j] , -1*lsp_y[1][i][j]
                    )
                );
            }
        }

        dt = CTS::CFL*min( dx/amax , dy/bmax );

        if( t+dt > time.second ){
            dt = time.second - t;
        }
    }

    // function to obtain the PLR of a given conserved variable
    // used as a subroutine in CU flux calculation
    vvvd get_plr(vvd& var, string var_name){
        int row = var.size();
        int col = var[0].size();
        
        vvvd plr; // [N, S, E, W] -> respectively

        vvd n(row, vd(col, 0));
        vvd s(row, vd(col, 0));
        vvd e(row, vd(col, 0));
        vvd w(row, vd(col, 0));

        for( int i=1 ; i<row-1 ; i++ ){
            for(int j=1 ; j<col-1 ; j++){
                double v1 = CTS::THETA * (var[i][j] - var[i][j-1]);
                double v2 = 0.5*(var[i][j+1] - var[i][j-1]);
                double v3 = CTS::THETA * (var[i][j+1] - var[i][j]);

                double slx = UTL::minmod( v1 , v2 , v3 );

                v1 = CTS::THETA * (var[i][j] - var[i-1][j]);
                v2 = 0.5*(var[i+1][j] - var[i-1][j]);
                v3 = CTS::THETA * (var[i+1][j] - var[i][j]);

                double sly = UTL::minmod( v1 , v2 , v3 );

                n[i][j] = var[i][j] + 0.5*sly;
                s[i][j] = var[i][j] - 0.5*sly;
                e[i][j] = var[i][j] + 0.5*slx;
                w[i][j] = var[i][j] - 0.5*slx;
            }
        }

        plr.push_back(n);
        plr.push_back(s);
        plr.push_back(e);
        plr.push_back(w);

        // PLR needs to be extended in a different manner than other extensions
        EXC::extend_plr(plr, bc, var_name);

        return plr;

    }

    // function to calculate the One sided local speeds of propagations
    // a+, a-
    vvvd get_lsp_a(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr){
        double GAMMA = CTS::GAMMA;

        // first we get the primitive variable values for East and West boundaries
        vvd uE = PRV::get_u(problem, rho_plr[2], mx_plr[2]);
        vvd uW = PRV::get_u(problem, rho_plr[3], mx_plr[3]);
        vvd vE = PRV::get_v(problem, rho_plr[2], my_plr[2]);
        vvd vW = PRV::get_v(problem, rho_plr[3], my_plr[3]);
        vvd pE = PRV::get_p(problem, rho_plr[2], uE, vE, E_plr[2]);        
        vvd pW = PRV::get_p(problem, rho_plr[3], uW, vW, E_plr[3]);        

        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();

        vvvd lsp;

        vvd ap(row, vd(col));
        vvd am(row, vd(col));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=0 ; j<col-1 ; j++){
                ap[i][j] = max(
                    0.0 ,
                    max(
                        uW[i][j+1] + sqrtf( abs((GAMMA*pW[i][j+1])/rho_plr[3][i][j+1]) ),
                        uE[i][j] + sqrtf( abs((GAMMA*pE[i][j]) / rho_plr[2][i][j]) )
                    )
                );

                am[i][j] = min(
                    0.0,
                    min(
                        uW[i][j+1] - sqrtf( abs((GAMMA*pW[i][j+1])/rho_plr[3][i][j+1]) ),
                        uE[i][j] - sqrtf( abs((GAMMA*pE[i][j])/rho_plr[2][i][j]) )
                    )
                );
            }
        }

        lsp.push_back(ap);
        lsp.push_back(am);

        return lsp;
    }   

    // b+, b-
    vvvd get_lsp_b(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr){
        double GAMMA = CTS::GAMMA;

        // first we get the primitive variable values for East and West boundaries
        vvd uN = PRV::get_u(problem, rho_plr[0], mx_plr[0]);
        vvd uS = PRV::get_u(problem, rho_plr[1], mx_plr[1]);
        vvd vN = PRV::get_v(problem, rho_plr[0], my_plr[0]);
        vvd vS = PRV::get_v(problem, rho_plr[1], my_plr[1]);
        vvd pN = PRV::get_p(problem, rho_plr[0], uN, vN, E_plr[0]);        
        vvd pS = PRV::get_p(problem, rho_plr[1], uS, vS, E_plr[1]);        

        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();

        vvvd lsp;

        vvd bp(row, vd(col));
        vvd bm(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                bp[i][j] = max(
                    0.0 ,
                    max(
                        vS[i+1][j] + sqrtf( abs((GAMMA*pS[i+1][j])/rho_plr[1][i+1][j]) ),
                        vN[i][j] + sqrtf( abs((GAMMA*pN[i][j]) / rho_plr[0][i][j]) )
                    )
                );

                bm[i][j] = min(
                    0.0,
                    min(
                        vS[i+1][j] - sqrtf( abs((GAMMA*pS[i+1][j])/rho_plr[1][i+1][j]) ),
                        vN[i][j] - sqrtf( abs((GAMMA*pN[i][j])/rho_plr[0][i][j]) )
                    )
                );
            }
        }

        lsp.push_back(bp);
        lsp.push_back(bm);

        return lsp;
    }

};

