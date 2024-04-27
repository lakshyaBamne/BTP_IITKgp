/*
    ! @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-8 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme - with RH optimization for Euler Equations of Gas Dynamics (CURH Scheme)
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>

#include "Constants.h"
#include "Utility.h"
#include "ExtendCells.h"
#include "OutputResult.h"
#include "PrimitiveVariables.h"

using namespace std;

namespace CTS = Constants;
namespace UTL = Utility;
namespace EXC = ExtendCells;
namespace OPR = OutputResult;
namespace PRV = PrimitiveVariables;

#define pdd pair<double,double>
#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define vvvvd vector< vector< vector< vector<double> > > >

class CURH{
public: // variables
    
    int problem; // [1-9]
    int mode; // [1-2]

    int Nx, Ny; // this is the number of intervals in the x and y directions
    pdd domainX, domainY;
    double dx, dy;
    pdd time; // [t_initial, t_final]
    double t, dt; // time at some instant of the simulations and time step

    vd X, Y; // X and Y grids => computational domain
    vvd rho, mx, my, E; // Conserved variable vectors

    unordered_map<string,string> bc;



public: // methods
    
    /*
        Class constructor is used to run the simulations based on the given variables
    */
    CURH(int problem, int mode) : problem(problem), mode(mode) {
        // first we need to initialize the variables used for the simulation
        init_vars();

        // now we can create the computational grids
        make_grids();

        // now we can initialize the conserved variables to start the simulations
        initialize_conserved_variables();

        // export the environment variables to an external file
        export_env(); 
    }

    /*
        Method to initialize the variables required to run the simulations
    */
    void init_vars();

    /*
        Method to initialize the computational grids based on the problem
    */
    void make_grids();

    /*
        Method to initialize the conserved variables to start the iterations
    */
    void initialize_conserved_variables();

    /*
        Method to export the environment variables related to the simulation to run
    */
    void export_env();

    /*
        ! Methods to run the CURH Scheme follow from here on
    */
    
    /*
        Method to run the CURH scheme and store the output after every 5 time steps
    */
    void RunCURH_complete();

    /*
        Method to run the CURH scheme and store the output for the initial 
        and final conditions only
    */
    void RunCURH_partial();

    /*
        !STEP-1 : Calculate the slopes of the Reconstructions (this is a second order scheme)
        !STEP-2 : Calculate the Piecewise Linear Reconstructions using the slopes 
    */
    vvvd get_plr(vvd& var, string var_name);

    vvd get_slx(vvd& var);
    vvd get_sly(vvd& var);

    /*
        !STEP-STAR : Method to update the conserved variable vectors at every time step
                    !using Three stage SSPRK time integration scheme
    */
    void update_conserved_variables();

    /*
        Method takes a set of conserved variable matrices and returns the CURH Numerical Flux 
        corresponding to it
    */
    vvvvd curh_flux(vvd rho, vvd mx, vvd my, vvd E, bool update_time);

    /*
        !STEP-3 : Calculate the CURH Numerical flux
    */
    // vvvd get_curh_x(vvvd& lspx);
    // vvvd get_curh_y(vvvd& lspy);

    vvvd get_curh_x(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr, vvvd& lspx);
    vvvd get_curh_y(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr, vvvd& lspy);

    /*
        Methods to get the Local speeds of propagation using the PLR
    */
    vvvd get_lsp_x(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr);
    vvvd get_lsp_y(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr);

    /*
        Method to update the time step using the CFL conditions
    */
    void update_dt(vvvd& lspx, vvvd& lspy);

};

/*
    Function Definitions
*/

void CURH::init_vars(){

    switch(this->mode){
        case 1: // CURH
            
            switch(this->problem){
                case 1: // MCW
                    
                    domainX.first = -0.2;
                    domainX.second = 0.2;

                    domainY.first = 0;
                    domainY.second = 0.8;

                    dx = (double)1/250;
                    dy = (double)1/250;

                    time.first = 0;
                    time.second = 2.0;

                    t = 0;

                    Nx = (int)( (domainX.second - domainX.first) / dx );
                    Ny = (int)( (domainY.second - domainY.first) / dy );

                    bc["N"] = "FREE";
                    bc["S"] = "FREE";
                    bc["E"] = "FREE";
                    bc["W"] = "FREE";

                    break;
                case 2: // 2DR-CFG3
                    
                    domainX.first = 0;
                    domainX.second = 1.2;

                    domainY.first = 0;
                    domainY.second = 1.2;

                    dx = (double)1/100;
                    dy = (double)1/100;

                    time.first = 0;
                    time.second = 1.0;

                    t = 0;

                    Nx = (int)( (domainX.second - domainX.first) / dx );
                    Ny = (int)( (domainY.second - domainY.first) / dy );

                    bc["N"] = "FREE";
                    bc["S"] = "FREE";
                    bc["E"] = "FREE";
                    bc["W"] = "FREE";

                    break;
                case 3: // EXP

                    domainX.first = 0;
                    domainX.second = 1.5;

                    domainY.first = 0;
                    domainY.second = 1.5;

                    dx = (double)3/800;
                    dy = (double)3/800;

                    time.first = 0;
                    time.second = 3.2;

                    t = 0;

                    Nx = (int)( (domainX.second - domainX.first) / dx );
                    Ny = (int)( (domainY.second - domainY.first) / dy );

                    bc["N"] = "FREE";
                    bc["S"] = "REF";
                    bc["E"] = "FREE";
                    bc["W"] = "REF";
                    
                    break;
                case 4: // IMP1
                    
                    break;
                case 5: // IMP2
                    
                    break;
                case 6: // KHI1
                    
                    break;
                case 7: // KHI2
                    
                    break;
                case 8: // KHI3
                    
                    break;
                case 9: // RTI
                    
                    break;
                
                default:
                    cout << "[ERROR] Please enter correct problem to run the simulations" << endl;
                    break;
            }

            break;
        case 2: // REFERENCE

            switch(this->problem){
                case 1: // MCW
                    
                    domainX.first = -0.2;
                    domainX.second = 0.2;

                    domainY.first = 0;
                    domainY.second = 0.8;

                    dx = (double)1/1000;
                    dy = (double)1/1000;

                    time.first = 0;
                    time.second = 2.0;

                    Nx = (int)( (domainX.second - domainX.first) / dx );
                    Ny = (int)( (domainY.second - domainY.first) / dy );

                    bc["N"] = "FREE";
                    bc["S"] = "FREE";
                    bc["E"] = "FREE";
                    bc["W"] = "FREE";

                    break;
                case 2: // 2DR-CFG3
                    
                    domainX.first = 0;
                    domainX.second = 1.2;

                    domainY.first = 0;
                    domainY.second = 1.2;

                    dx = (double)1/4000;
                    dy = (double)1/4000;

                    time.first = 0;
                    time.second = 1.0;

                    Nx = (int)( (domainX.second - domainX.first) / dx );
                    Ny = (int)( (domainY.second - domainY.first) / dy );

                    bc["N"] = "FREE";
                    bc["S"] = "FREE";
                    bc["E"] = "FREE";
                    bc["W"] = "FREE";

                    break;
                case 3: // EXP

                    domainX.first = 0;
                    domainX.second = 1.5;

                    domainY.first = 0;
                    domainY.second = 1.5;

                    dx = (double)3/3200;
                    dy = (double)3/3200;

                    time.first = 0;
                    time.second = 3.2;

                    Nx = (int)( (domainX.second - domainX.first) / dx );
                    Ny = (int)( (domainY.second - domainY.first) / dy );

                    bc["N"] = "FREE";
                    bc["S"] = "REF";
                    bc["E"] = "FREE";
                    bc["W"] = "REF";
                    
                    break;
                case 4: // IMP1
                    
                    break;
                case 5: // IMP2
                    
                    break;
                case 6: // KHI1
                    
                    break;
                case 7: // KHI2
                    
                    break;
                case 8: // KHI3
                    
                    break;
                case 9: // RTI
                    
                    break;
                
                default:
                    cout << "[ERROR] Please enter correct problem to run the simulations" << endl;
                    break;
            }

            break;

        default:
            cout << "[ERROR] Please enter correct mode of operation" << endl;
            break;
    }

    // log message
    cout << "---Initialized Simulation Variables---" << endl;

}

void CURH::make_grids(){

    // X
    for(int i=0 ; i<Nx ; i++){
        X.push_back(domainX.first + (0.5+i)*dx);
    }

    // Y
    for(int i=0 ; i<Ny ; i++){
        Y.push_back(domainY.first + (0.5+i)*dy);
    }

    // Log message
    cout << "---Initialized Computational Grid---" << endl;
    cout << "X : ";
    for(int i=0 ; i<10 ; i++){
        cout << X[i] << ", ";
    }
    cout << "..., " << X[X.size()-2] << ", " << X.back();
    cout << endl;

    cout << "Y : ";
    for(int i=0 ; i<10 ; i++){
        cout << Y[i] << ", ";
    }
    cout << "..., " << Y[Y.size()-2] << ", " << Y.back();
    cout << endl;

    // write the computational grids to an external file
    OPR::write_grids(mode, X, Y);
}

void CURH::initialize_conserved_variables(){

    switch(this->problem){
        case 1: // MCW
            
            for(auto j : Y){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : X){
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

            break;
        case 2: // 2DR-CFG3

            for(auto j : Y){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : X){
                    
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
            
            break;
        case 3: // EXP

            for(auto j : Y){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : X){
                    
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
            
            break;
        case 4: // IMP1
            
            break;
        case 5: // IMP2
            
            break;
        case 6: // KHI1
            
            break;
        case 7: // KHI2
        
            break;
        case 8: // KHI3
            
            break;
        case 9: // RTI
            
            break;
        
        default:
            cout << "[ERROR] Please enter the correct problem" << endl;
            break;
    }

    cout << "---Initialized Conserved variables---" << endl;

    // Now we need to extend the conserved variables according to the boundary conditions
    EXC::extend_matrix(rho, bc, "Density");
    EXC::extend_matrix(mx, bc, "MomentumX");
    EXC::extend_matrix(my, bc, "MomentumY");
    EXC::extend_matrix(E, bc, "Energy");

    cout << "---Extended cells for the conserved variables---" << endl;


}

void CURH::export_env(){
    string file_name;

    if( mode == 1 ){
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
    ss << domainX.first << " " << domainX.second << endl;
    ss << domainY.first << " " << domainY.second << endl;
    ss << time.first << " " << time.second << endl;
    ss << Nx << endl;
    ss << Ny << endl;

    fout << ss.str() << endl;

    fout.close();
}

/*
    !Methods to run the CURH Scheme after the environment has been made
*/

void CURH::RunCURH_complete(){

    vvd u;
    vvd v;
    vvd p;

    while( t < time.second ){
        cout << "t = " << t << " | dt = " << dt << endl;
        
        //! write variable values at each step for making animations
        OPR::write_matrix(mode, "Density", rho);
        OPR::write_matrix(mode, "MomentumX", mx);
        OPR::write_matrix(mode, "MomentumY", my);
        OPR::write_matrix(mode, "Energy", E);

        // get the primitive variables
        u = PRV::get_u(problem, rho, mx);
        v = PRV::get_v(problem, rho, my);
        p = PRV::get_p(problem, rho, u, v, E);

        OPR::write_matrix(mode, "VelocityX", u);
        OPR::write_matrix(mode, "VelocityY", v);
        OPR::write_matrix(mode, "Pressure", p);

        // function to update the conservative variables internally
        // also update dt for next iteration
        update_conserved_variables();

        // update new time step calculated using CFL conditions
        t += dt;

    }
}

void CURH::RunCURH_partial(){
    //! write initial values for plotting
    OPR::write_matrix(mode, "Density", rho);
    OPR::write_matrix(mode, "MomentumX", mx);
    OPR::write_matrix(mode, "MomentumY", my);
    OPR::write_matrix(mode, "Energy", E);  

    while( t < time.second ){
        cout << "t = " << t << " | dt = " << dt << endl;
        
        // function to update the conservative variables internally
        // also update dt for next iteration
        update_conserved_variables();

        // update new time step calculated using CFL conditions
        t += dt;

    }

    //! write final values for plotting
    OPR::write_matrix(mode, "Density", rho);
    OPR::write_matrix(mode, "MomentumX", mx);
    OPR::write_matrix(mode, "MomentumY", my);
    OPR::write_matrix(mode, "Energy", E);
}

vvvd CURH::get_plr(vvd& var, string var_name){
    vvd slx = get_slx(var);        
    vvd sly = get_sly(var);        

    int row = var.size();
    int col = var[0].size();
    
    vvvd plr; // [N, S, E, W, NE, NW, SE, SW] -> respectively

    vvd n(row, vd(col, 0));
    vvd s(row, vd(col, 0));
    vvd e(row, vd(col, 0));
    vvd w(row, vd(col, 0));
    vvd ne(row, vd(col, 0));
    vvd nw(row, vd(col, 0));
    vvd se(row, vd(col, 0));
    vvd sw(row, vd(col, 0));

    for( int i=1 ; i<row-1 ; i++ ){
        for(int j=1 ; j<col-1 ; j++){
            n[i][j] = var[i][j] + 0.5*sly[i][j];
            s[i][j] = var[i][j] - 0.5*sly[i][j];
            e[i][j] = var[i][j] + 0.5*slx[i][j];
            w[i][j] = var[i][j] - 0.5*slx[i][j];
            ne[i][j] = var[i][j] + 0.5*slx[i][j] + 0.5*sly[i][j];
            nw[i][j] = var[i][j] - 0.5*slx[i][j] + 0.5*sly[i][j];
            se[i][j] = var[i][j] + 0.5*slx[i][j] - 0.5*sly[i][j];
            sw[i][j] = var[i][j] - 0.5*slx[i][j] - 0.5*sly[i][j];
        }
    }

    plr.push_back(n);
    plr.push_back(s);
    plr.push_back(e);
    plr.push_back(w);
    plr.push_back(ne);
    plr.push_back(nw);
    plr.push_back(se);
    plr.push_back(sw);

    // extend cells for the piecewise linear reconstructions
    EXC::extend_plr(plr, bc, var_name);

    return plr;

}

vvd CURH::get_slx(vvd& var){
    int row = var.size();
    int col = var[0].size();

    vvd slx(row, vd(col));

    for(int i=1 ; i<row-1 ; i++){
        for(int j=1 ; j<col-1 ; j++){
            double v1 = CTS::THETA * (var[i][j] - var[i][j-1]);
            double v2 = 0.5*(var[i][j+1] - var[i][j-1]);
            double v3 = CTS::THETA * (var[i][j+1] - var[i][j]);

            slx[i][j] = UTL::minmod( v1 , v2 , v3 );
        }
    }

    return slx;
}

vvd CURH::get_sly(vvd& var){
    int row = var.size();
    int col = var[0].size();

    vvd sly(row, vd(col));

    for(int i=1 ; i<row-1 ; i++){
        for(int j=1 ; j<col-1 ; j++){
            double v1 = CTS::THETA * (var[i][j] - var[i-1][j]);
            double v2 = 0.5*(var[i+1][j] - var[i-1][j]);
            double v3 = CTS::THETA * (var[i+1][j] - var[i][j]);

            sly[i][j] = UTL::minmod( v1 , v2 , v3 );
        }
    }

    return sly;
}

// vvvd CURH::get_curh_x(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr, vvvd& lspx){

//     int row = rho_plr[0].size();
//     int col = rho_plr[0][0].size();

//     // CURH Numerical flux
//     vvd F1(row, vd(col,0));
//     vvd F2(row, vd(col,0));
//     vvd F3(row, vd(col,0));
//     vvd F4(row, vd(col,0));

//     for(int i=1 ; i<row-1 ; i++){
//         for(int j=0 ; j<col-1 ; j++){
//             // Primitive Variables
//             double uE = mx_plr[2][i][j] / rho_plr[2][i][j];
//             double uW = mx_plr[3][i][j+1] / rho_plr[3][i][j+1];
//             double vE = my_plr[2][i][j] / rho_plr[2][i][j];
//             double vW = my_plr[3][i][j+1] / rho_plr[3][i][j+1];
//             double pE = ( CTS::GAMMA - 1 ) * ( E_plr[2][i][j] - 0.5*rho_plr[2][i][j]*( pow(uE,2) + pow(vE,2) ) );
//             double pW = ( CTS::GAMMA - 1 ) * ( E_plr[3][i][j+1] - 0.5*rho_plr[3][i][j+1]*( pow(uW,2) + pow(vW,2) ) );        

//             // Flux vector elements
//             double f1E = mx_plr[2][i][j];
//             double f2E = mx_plr[2][i][j]*uE + pE;
//             double f3E = mx_plr[2][i][j]*vE;
//             double f4E = uE*( E_plr[2][i][j] + pE );

//             double f1W = mx_plr[3][i][j+1];
//             double f2W = mx_plr[3][i][j+1]*uW + pW;
//             double f3W = vW*mx_plr[3][i][j+1];
//             double f4W = uW*( E_plr[3][i][j+1] + pW );

//             double dist = lspx[0][i][j] - lspx[1][i][j];
//             double prod = lspx[0][i][j]*lspx[1][i][j];

//             if( dist > CTS::EPSILON ){
//                 // intermediate terms
//                 double rho_star = ( lspx[0][i][j]*rho_plr[3][i][j+1] - lspx[1][i][j]*rho_plr[2][i][j] - (f1W-f1E) ) / dist;
//                 double mx_star = ( lspx[0][i][j]*mx_plr[3][i][j+1] - lspx[1][i][j]*mx_plr[2][i][j] - (f2W-f2E) ) / dist;
//                 double my_star = ( lspx[0][i][j]*my_plr[3][i][j+1] - lspx[1][i][j]*my_plr[2][i][j] - (f3W-f3E) ) / dist;
//                 double E_star = ( lspx[0][i][j]*E_plr[3][i][j+1] - lspx[1][i][j]*E_plr[2][i][j] - (f4W-f4E) ) / dist;

//                 double drho = UTL::minmod(
//                     rho_star - rho_plr[4][i][j],
//                     rho_plr[5][i][j+1] - rho_star,
//                     rho_star-rho_plr[6][i][j],
//                     rho_plr[7][i][j+1]-rho_star
//                 );

//                 double dmx = UTL::minmod(
//                     mx_star - mx_plr[4][i][j],
//                     mx_plr[5][i][j+1] - mx_star,
//                     mx_star-mx_plr[6][i][j],
//                     mx_plr[7][i][j+1]-mx_star
//                 );

//                 double dmy = UTL::minmod(
//                     my_star - my_plr[4][i][j],
//                     my_plr[5][i][j+1] - my_star,
//                     my_star-my_plr[6][i][j],
//                     my_plr[7][i][j+1]-my_star
//                 );

//                 double dE = UTL::minmod(
//                     E_star - E_plr[4][i][j],
//                     E_plr[5][i][j+1] - E_star,
//                     E_star-E_plr[6][i][j],
//                     E_plr[7][i][j+1]-E_star
//                 );

//                 F1[i][j] = ( (lspx[0][i][j]*f1E - lspx[1][i][j]*f1W) + prod*( rho_plr[3][i][j+1] - rho_plr[2][i][j] - drho ) ) / dist;
//                 F2[i][j] = ( (lspx[0][i][j]*f2E - lspx[1][i][j]*f2W) + prod*( mx_plr[3][i][j+1] - mx_plr[2][i][j] - dmx ) ) / dist;
//                 F3[i][j] = ( (lspx[0][i][j]*f3E - lspx[1][i][j]*f3W) + prod*( my_plr[3][i][j+1] - my_plr[2][i][j] - dmy ) ) / dist;
//                 F4[i][j] = ( (lspx[0][i][j]*f4E - lspx[1][i][j]*f4W) + prod*( E_plr[3][i][j+1] - E_plr[2][i][j] - dE ) ) / dist;
//             }
//             else{
//                 F1[i][j] = 0.5*( f1E + f1W );
//                 F2[i][j] = 0.5*( f2E + f2W );
//                 F3[i][j] = 0.5*( f3E + f3W );
//                 F4[i][j] = 0.5*( f4E + f4W );
//             }
//         }
//     }

//     return {F1, F2, F3, F4};

// }

// vvvd CURH::get_curh_y(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr, vvvd& lspy){

//     int row = rho_plr[0].size();
//     int col = rho_plr[0][0].size();

//     // CURH Numerical flux
//     vvd G1(row, vd(col,0));
//     vvd G2(row, vd(col,0));
//     vvd G3(row, vd(col,0));
//     vvd G4(row, vd(col,0));

//     for(int i=0 ; i<row-1 ; i++){
//         for(int j=1 ; j<col-1 ; j++){
//             // first we need to calculate the required primitive variables
//             double uN = mx_plr[0][i][j] / rho_plr[0][i][j];
//             double uS = mx_plr[1][i+1][j] / rho_plr[1][i+1][j];
//             double vN = my_plr[0][i][j] / rho_plr[0][i][j];
//             double vS = my_plr[1][i+1][j] / rho_plr[1][i+1][j];
//             double pN = ( CTS::GAMMA - 1 )*( E_plr[0][i][j] - 0.5*rho_plr[0][i][j]*( pow(uN, 2) + pow(vN,2) ) );
//             double pS = ( CTS::GAMMA - 1 )*( E_plr[1][i+1][j] - 0.5*rho_plr[1][i+1][j]*( pow(uS, 2) + pow(vS,2) ) );

//             double g1N = my_plr[0][i][j];
//             double g2N = my_plr[0][i][j]*uN;
//             double g3N = my_plr[0][i][j]*vN + pN;
//             double g4N = vN*( E_plr[0][i][j] + pN );

//             double g1S = my_plr[1][i+1][j];
//             double g2S = my_plr[1][i+1][j]*uS;
//             double g3S = my_plr[1][i+1][j]*vS + pS;
//             double g4S = vS*( E_plr[1][i+1][j] + pS );

//             double dist = lspy[0][i][j] - lspy[1][i][j];
//             double prod = lspy[0][i][j]*lspy[1][i][j];

//             if( dist > CTS::EPSILON ){
//                 double rho_star = ( lspy[0][i][j]*rho_plr[1][i+1][j] - lspy[1][i][j]*rho_plr[0][i][j] - (g1S-g1N) ) / dist;
//                 double mx_star = ( lspy[0][i][j]*mx_plr[1][i+1][j] - lspy[1][i][j]*mx_plr[0][i][j] - (g2S-g2N) ) / dist;
//                 double my_star = ( lspy[0][i][j]*my_plr[1][i+1][j] - lspy[1][i][j]*my_plr[0][i][j] - (g3S-g3N) ) / dist;
//                 double E_star = ( lspy[0][i][j]*E_plr[1][i+1][j] - lspy[1][i][j]*E_plr[0][i][j] - (g4S-g4N) ) / dist;

//                 double drho = UTL::minmod(
//                     rho_star - rho_plr[4][i][j],
//                     rho_plr[6][i+1][j] - rho_star,
//                     rho_star - rho_plr[5][i][j],
//                     rho_plr[7][i+1][j] - rho_star
//                 );

//                 double dmx = UTL::minmod(
//                     mx_star - mx_plr[4][i][j],
//                     mx_plr[6][i+1][j] - mx_star,
//                     mx_star - mx_plr[5][i][j],
//                     mx_plr[7][i+1][j] - mx_star
//                 );

//                 double dmy = UTL::minmod(
//                     my_star - my_plr[4][i][j],
//                     my_plr[6][i+1][j] - my_star,
//                     my_star - my_plr[5][i][j],
//                     my_plr[7][i+1][j] - my_star
//                 );

//                 double dE = UTL::minmod(
//                     E_star - E_plr[4][i][j],
//                     E_plr[6][i+1][j] - E_star,
//                     E_star - E_plr[5][i][j],
//                     E_plr[7][i+1][j] - E_star
//                 );

//                 G1[i][j] = ( (lspy[0][i][j]*g1N - lspy[1][i][j]*g1S) + prod*( ( rho_plr[1][i+1][j] - rho_plr[0][i][j] ) - drho ) ) / dist;
//                 G2[i][j] = ( (lspy[0][i][j]*g2N - lspy[1][i][j]*g2S) + prod*( ( mx_plr[1][i+1][j] - mx_plr[0][i][j] ) - dmx ) ) / dist;
//                 G3[i][j] = ( (lspy[0][i][j]*g3N - lspy[1][i][j]*g3S) + prod*( ( my_plr[1][i+1][j] - my_plr[0][i][j] ) - dmy ) ) / dist;
//                 G4[i][j] = ( (lspy[0][i][j]*g4N - lspy[1][i][j]*g4S) + prod*( ( E_plr[1][i+1][j] - E_plr[0][i][j] ) - dE ) ) / dist;
//             }  
//             else{
//                 G1[i][j] = 0.5*( g1N + g1S );                
//                 G2[i][j] = 0.5*( g2N + g2S );                
//                 G3[i][j] = 0.5*( g3N + g3S );                
//                 G4[i][j] = 0.5*( g4N + g4S );                
//             }

//         }
//     }

//     return {G1, G2, G3, G4};
// }

// vvvd CURH::get_lsp_x(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr){
//     // first coordinate represents the row and the y-index
//     int row = rho_plr[0].size();
//     int col = rho_plr[0][0].size();

//     vvd ap(row, vd(col,0));
//     vvd am(row, vd(col,0));

//     for(int i=1 ; i<row-1 ; i++){
//         for(int j=0 ; j<col-1 ; j++){
//             // first we need to calculate the required primitive variables
//             double uE = mx_plr[2][i][j] / rho_plr[2][i][j];
//             double uW = mx_plr[3][i][j+1] / rho_plr[3][i][j+1];
//             double vE = my_plr[2][i][j] / rho_plr[2][i][j];
//             double vW = my_plr[3][i][j+1] / rho_plr[3][i][j+1];
//             double pE = ( CTS::GAMMA - 1 ) * ( E_plr[2][i][j] - 0.5*rho_plr[2][i][j]*( pow(uE,2) + pow(vE,2) ) );
//             double pW = ( CTS::GAMMA - 1 ) * ( E_plr[3][i][j+1] - 0.5*rho_plr[3][i][j+1]*( pow(uW,2) + pow(vW,2) ) );

//             // now we can calulcate the Local Speeds of propagation
//             ap[i][j] = max(
//                 0.0,
//                 max(
//                     uW + sqrt( CTS::GAMMA*pW / rho_plr[3][i][j+1] ),
//                     uE + sqrt( CTS::GAMMA*pE / rho_plr[2][i][j] )
//                 )
//             );

//             am[i][j] = min(
//                 0.0,
//                 min(
//                     uW - sqrt( CTS::GAMMA*pW/rho_plr[3][i][j+1] ),
//                     uE - sqrt( CTS::GAMMA*pE/rho_plr[2][i][j] )
//                 )
//             );

//             //! The following optimization on the local speeds of propagation is part of CURH Scheme
//             //! if it is excluded then the code runs for CU scheme
//             if( CTS::ScaleSw ){
//                 double differ1 = abs(0.5*( rho_plr[3][i][j+1]*pow(uW,2) - rho_plr[2][i][j]*pow(uE,2) ) + ( (pW-pE) / (CTS::GAMMA-1) ));
//                 double differ2 = 0.5*abs( rho_plr[3][i][j+1]*pow(vW,2) - rho_plr[2][i][j]*pow(vE,2) );
//                 double differ = sqrt(pow(differ1,2)+pow(differ2,2));

//                 double alpha;

//                 if( differ > CTS::EPSILON ){
//                     alpha = differ1/differ;
//                 }
//                 else{
//                     alpha = 0;
//                 }

//                 ap[i][j] = max(
//                     0.0,
//                     max(
//                         uW + alpha*sqrt(CTS::GAMMA*pW/rho_plr[3][i][j+1]),
//                         uE + alpha*sqrt(CTS::GAMMA*pE/rho_plr[2][i][j])
//                     )
//                 );

//                 am[i][j] = min(
//                     0.0,
//                     min(
//                         uW - alpha*sqrt(CTS::GAMMA*pW/rho_plr[3][i][j+1]),
//                         uE - alpha*sqrt(CTS::GAMMA*pE/rho_plr[2][i][j])
//                     )
//                 );
//             }

//             double f1E = mx_plr[2][i][j];
//             double f2E = mx_plr[2][i][j]*uE + pE;
//             double f3E = mx_plr[2][i][j]*vE;
//             double f4E = uE*( E_plr[2][i][j] + pE );

//             double f1W = mx_plr[3][i][j+1];
//             double f2W = mx_plr[3][i][j+1]*uW + pW;
//             double f3W = vW*mx_plr[3][i][j+1];
//             double f4W = uW*( E_plr[3][i][j+1] + pW );

//             if( CTS::ShrinkSw ){
//                 double df1 = f1W - f1E;
//                 double df2 = f2W - f2E;
//                 double df3 = f3W - f3E;
//                 double df4 = f4W - f4E;

//                 double du1 = rho_plr[3][i][j+1] - rho_plr[2][i][j];
//                 double du2 = mx_plr[3][i][j+1] - mx_plr[2][i][j];
//                 double du3 = my_plr[3][i][j+1] - my_plr[2][i][j];
//                 double du4 = E_plr[3][i][j+1] - E_plr[2][i][j];

//                 double du1eps, du2eps, du3eps, du4eps;

//                 if( du1 > 0 ){
//                     du1eps = max(du1, 1.0E-10);
//                 }
//                 else{  
//                     du1eps = min(du1, -1.0E-10);
//                 }

//                 if( du2 > 0 ){
//                     du2eps = max(du2, 1.0E-10);
//                 }
//                 else{  
//                     du2eps = min(du2, -1.0E-10);
//                 }

//                 if( du3 > 0 ){
//                     du3eps = max(du3, 1.0E-10);
//                 }
//                 else{  
//                     du3eps = min(du3, -1.0E-10);
//                 }

//                 if( du4 > 0 ){
//                     du4eps = max(du4, 1.0E-10);
//                 }
//                 else{  
//                     du4eps = min(du4, -1.0E-10);
//                 }

//                 double ratio1 = 2.0*df1 / (du1+du1eps);
//                 double ratio2 = 2.0*df2 / (du2+du2eps);
//                 double ratio3 = 2.0*df3 / (du3+du3eps);
//                 double ratio4 = 2.0*df4 / (du4+du4eps);

//                 double ratiomin = min(ratio1, min(ratio2, min(ratio3, ratio4)));
//                 double ratiomax = max(ratio1, max(ratio2, max(ratio3, ratio4)));

//                 if( ratiomax > 0 ){
//                     ap[i][j] = min( ap[i][j] , ratiomax );
//                     am[i][j] = max( am[i][j] , -ratiomax );
//                 }

//                 if( ratiomin < 0 ){
//                     ap[i][j] = min( ap[i][j] , -ratiomin );
//                     am[i][j] = max( am[i][j] , ratiomin );
//                 }

//             }            

//         }
//     }

//     return {ap,am};
// }

// vvvd CURH::get_lsp_y(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr){
//     // first coordinate represents the row and the y-index
//     int row = rho_plr[0].size();
//     int col = rho_plr[0][0].size();

//     vvd bp(row, vd(col,0));
//     vvd bm(row, vd(col,0));

//     for(int i=0 ; i<row-1 ; i++){
//         for(int j=1 ; j<col-1 ; j++){
//             // first we need to calculate the required primitive variables
//             double uN = mx_plr[0][i][j] / rho_plr[0][i][j];
//             double uS = mx_plr[1][i+1][j] / rho_plr[1][i+1][j];
//             double vN = my_plr[0][i][j] / rho_plr[0][i][j];
//             double vS = my_plr[1][i+1][j] / rho_plr[1][i+1][j];
//             double pN = ( CTS::GAMMA - 1 )*( E_plr[0][i][j] - 0.5*rho_plr[0][i][j]*( pow(uN, 2) + pow(vN,2) ) );
//             double pS = ( CTS::GAMMA - 1 )*( E_plr[1][i+1][j] - 0.5*rho_plr[1][i+1][j]*( pow(uS, 2) + pow(vS,2) ) );
        
//             // now we can calulcate the Local Speeds of propagation
//             bp[i][j] = max(
//                 0.0,
//                 max(
//                     vS + sqrt( CTS::GAMMA*pS/rho_plr[1][i+1][j] ),
//                     vN + sqrt( CTS::GAMMA*pN/rho_plr[0][i][j] )
//                 )
//             );

//             bm[i][j] = min(
//                 0.0,
//                 min(
//                     vS - sqrt( CTS::GAMMA*pS/rho_plr[1][i+1][j] ),
//                     vN - sqrt( CTS::GAMMA*pN/rho_plr[0][i][j] )
//                 )
//             );

//             //! The following optimization on the local speeds of propagation is part of CURH Scheme
//             //! if it is excluded then the code runs for CU scheme
//             if( CTS::ScaleSw ){
//                 double differ1 = abs(0.5*( rho_plr[1][i+1][j]*pow(vS,2) - rho_plr[0][i][j]*pow(vN,2) ) + ( (pS-pN) / (CTS::GAMMA-1) ));
//                 double differ2 = 0.5*abs( rho_plr[1][i+1][j]*pow(uS,2) - rho_plr[0][i][j]*pow(uN,2) );
//                 double differ = sqrt(pow(differ1,2)+pow(differ2,2));

//                 double alpha;

//                 if( differ > CTS::EPSILON ){
//                     alpha = differ1/differ;
//                 }
//                 else{
//                     alpha = 0;
//                 }

//                 bp[i][j] = max(
//                     0.0,
//                     max(
//                         vS + alpha*sqrt(CTS::GAMMA*pS/rho_plr[1][i+1][j]),
//                         vN + alpha*sqrt(CTS::GAMMA*pN/rho_plr[0][i][j])
//                     )
//                 );

//                 bm[i][j] = min(
//                     0.0,
//                     min(
//                         vS - alpha*sqrt(CTS::GAMMA*pS/rho_plr[1][i+1][j]),
//                         vN - alpha*sqrt(CTS::GAMMA*pN/rho_plr[0][i][j])
//                     )
//                 );
//             }

//             double g1N = my_plr[0][i][j];
//             double g2N = my_plr[0][i][j]*uN;
//             double g3N = my_plr[0][i][j]*vN + pN;
//             double g4N = vN*( E_plr[0][i][j] + pN );

//             double g1S = my_plr[1][i+1][j];
//             double g2S = my_plr[1][i+1][j]*uS;
//             double g3S = my_plr[1][i+1][j]*vS + pS;
//             double g4S = vS*( E_plr[1][i+1][j] + pS );

//             if( CTS::ShrinkSw ){
//                 double dg1 = g1S - g1N;
//                 double dg2 = g2S - g2N;
//                 double dg3 = g3S - g3N;
//                 double dg4 = g4S - g4N;

//                 double du1 = rho_plr[1][i+1][j] - rho_plr[0][i][j];
//                 double du2 = mx_plr[1][i+1][j] - mx_plr[0][i][j];
//                 double du3 = my_plr[1][i+1][j] - my_plr[0][i][j];
//                 double du4 = E_plr[1][i+1][j] - E_plr[0][i][j];

//                 double du1eps, du2eps, du3eps, du4eps;

//                 if( du1 > 0 ){
//                     du1eps = max(du1, 1.0E-10);
//                 }
//                 else{  
//                     du1eps = min(du1, -1.0E-10);
//                 }

//                 if( du2 > 0 ){
//                     du2eps = max(du2, 1.0E-10);
//                 }
//                 else{  
//                     du2eps = min(du2, -1.0E-10);
//                 }

//                 if( du3 > 0 ){
//                     du3eps = max(du3, 1.0E-10);
//                 }
//                 else{  
//                     du3eps = min(du3, -1.0E-10);
//                 }

//                 if( du4 > 0 ){
//                     du4eps = max(du4, 1.0E-10);
//                 }
//                 else{  
//                     du4eps = min(du4, -1.0E-10);
//                 }

//                 double ratio1 = 2.0*dg1 / (du1+du1eps);
//                 double ratio2 = 2.0*dg2 / (du2+du2eps);
//                 double ratio3 = 2.0*dg3 / (du3+du3eps);
//                 double ratio4 = 2.0*dg4 / (du4+du4eps);

//                 double ratiomin = min(ratio1, min(ratio2, min(ratio3, ratio4)));
//                 double ratiomax = max(ratio1, max(ratio2, max(ratio3, ratio4)));

//                 if( ratiomax > 0 ){
//                     bp[i][j] = min( bp[i][j] , ratiomax );
//                     bm[i][j] = max( bm[i][j] , -ratiomax );
//                 }

//                 if( ratiomin < 0 ){
//                     bp[i][j] = min( bp[i][j] , -ratiomin );
//                     bm[i][j] = max( bm[i][j] , ratiomin );
//                 }
//             }

//         }
//     }

//     return {bp,bm};

// }

vvvvd CURH::curh_flux(vvd rho, vvd mx, vvd my, vvd E, bool update_time){
    // Calculate the Piecewise Linear Reconstruction
    vvvd rho_plr = get_plr(rho, "Density");
    vvvd mx_plr = get_plr(mx, "MomentumX");
    vvvd my_plr = get_plr(my, "MomentumY");
    vvvd E_plr = get_plr(E, "Energy");    

    int row = rho_plr[0].size();
    int col = rho_plr[0][0].size();

    //! lspx calculation
    vvd ap(row, vd(col,0));
    vvd am(row, vd(col,0));

    // CURH Numerical flux
    vvd F1(row, vd(col,0));
    vvd F2(row, vd(col,0));
    vvd F3(row, vd(col,0));
    vvd F4(row, vd(col,0));

    for(int i=1 ; i<row-1 ; i++){
        for(int j=0 ; j<col-1 ; j++){
            // first we need to calculate the required primitive variables
            double uE = mx_plr[2][i][j] / rho_plr[2][i][j];
            double uW = mx_plr[3][i][j+1] / rho_plr[3][i][j+1];
            double vE = my_plr[2][i][j] / rho_plr[2][i][j];
            double vW = my_plr[3][i][j+1] / rho_plr[3][i][j+1];
            double pE = ( CTS::GAMMA - 1 ) * ( E_plr[2][i][j] - 0.5*rho_plr[2][i][j]*( pow(uE,2) + pow(vE,2) ) );
            double pW = ( CTS::GAMMA - 1 ) * ( E_plr[3][i][j+1] - 0.5*rho_plr[3][i][j+1]*( pow(uW,2) + pow(vW,2) ) );

            // now we can calulcate the Local Speeds of propagation
            ap[i][j] = max(
                0.0,
                max(
                    uW + sqrt( CTS::GAMMA*pW / rho_plr[3][i][j+1] ),
                    uE + sqrt( CTS::GAMMA*pE / rho_plr[2][i][j] )
                )
            );

            am[i][j] = min(
                0.0,
                min(
                    uW - sqrt( CTS::GAMMA*pW/rho_plr[3][i][j+1] ),
                    uE - sqrt( CTS::GAMMA*pE/rho_plr[2][i][j] )
                )
            );

            //! The following optimization on the local speeds of propagation is part of CURH Scheme
            //! if it is excluded then the code runs for CU scheme
            if( CTS::ScaleSw ){
                double differ1 = abs(0.5*( rho_plr[3][i][j+1]*pow(uW,2) - rho_plr[2][i][j]*pow(uE,2) ) + ( (pW-pE) / (CTS::GAMMA-1) ));
                double differ2 = 0.5*abs( rho_plr[3][i][j+1]*pow(vW,2) - rho_plr[2][i][j]*pow(vE,2) );
                double differ = sqrt(pow(differ1,2)+pow(differ2,2));

                double alpha;

                if( differ > CTS::EPSILON ){
                    alpha = differ1/differ;
                }
                else{
                    alpha = 0;
                }

                ap[i][j] = max(
                    0.0,
                    max(
                        uW + alpha*sqrt(CTS::GAMMA*pW/rho_plr[3][i][j+1]),
                        uE + alpha*sqrt(CTS::GAMMA*pE/rho_plr[2][i][j])
                    )
                );

                am[i][j] = min(
                    0.0,
                    min(
                        uW - alpha*sqrt(CTS::GAMMA*pW/rho_plr[3][i][j+1]),
                        uE - alpha*sqrt(CTS::GAMMA*pE/rho_plr[2][i][j])
                    )
                );
            }

            double f1E = mx_plr[2][i][j];
            double f2E = mx_plr[2][i][j]*uE + pE;
            double f3E = mx_plr[2][i][j]*vE;
            double f4E = uE*( E_plr[2][i][j] + pE );

            double f1W = mx_plr[3][i][j+1];
            double f2W = mx_plr[3][i][j+1]*uW + pW;
            double f3W = vW*mx_plr[3][i][j+1];
            double f4W = uW*( E_plr[3][i][j+1] + pW );

            if( CTS::ShrinkSw ){
                double df1 = f1W - f1E;
                double df2 = f2W - f2E;
                double df3 = f3W - f3E;
                double df4 = f4W - f4E;

                double du1 = rho_plr[3][i][j+1] - rho_plr[2][i][j];
                double du2 = mx_plr[3][i][j+1] - mx_plr[2][i][j];
                double du3 = my_plr[3][i][j+1] - my_plr[2][i][j];
                double du4 = E_plr[3][i][j+1] - E_plr[2][i][j];

                double du1eps, du2eps, du3eps, du4eps;

                if( du1 > 0 ){
                    du1eps = max(du1, 1.0E-10);
                }
                else{  
                    du1eps = min(du1, -1.0E-10);
                }

                if( du2 > 0 ){
                    du2eps = max(du2, 1.0E-10);
                }
                else{  
                    du2eps = min(du2, -1.0E-10);
                }

                if( du3 > 0 ){
                    du3eps = max(du3, 1.0E-10);
                }
                else{  
                    du3eps = min(du3, -1.0E-10);
                }

                if( du4 > 0 ){
                    du4eps = max(du4, 1.0E-10);
                }
                else{  
                    du4eps = min(du4, -1.0E-10);
                }

                double ratio1 = 2.0*df1 / (du1+du1eps);
                double ratio2 = 2.0*df2 / (du2+du2eps);
                double ratio3 = 2.0*df3 / (du3+du3eps);
                double ratio4 = 2.0*df4 / (du4+du4eps);

                double ratiomin = min(ratio1, min(ratio2, min(ratio3, ratio4)));
                double ratiomax = max(ratio1, max(ratio2, max(ratio3, ratio4)));

                if( ratiomax > 0 ){
                    ap[i][j] = min( ap[i][j] , ratiomax );
                    am[i][j] = max( am[i][j] , -ratiomax );
                }

                if( ratiomin < 0 ){
                    ap[i][j] = min( ap[i][j] , -ratiomin );
                    am[i][j] = max( am[i][j] , ratiomin );
                }

            }         

            // calculation of CURH Flux in x direction starts
            double dist = ap[i][j] - am[i][j];
            double prod = ap[i][j]*am[i][j];

            if( dist > CTS::EPSILON ){
                // intermediate terms
                double rho_star = ( ap[i][j]*rho_plr[3][i][j+1] - am[i][j]*rho_plr[2][i][j] - (f1W-f1E) ) / dist;
                double mx_star = ( ap[i][j]*mx_plr[3][i][j+1] - am[i][j]*mx_plr[2][i][j] - (f2W-f2E) ) / dist;
                double my_star = ( ap[i][j]*my_plr[3][i][j+1] - am[i][j]*my_plr[2][i][j] - (f3W-f3E) ) / dist;
                double E_star = ( ap[i][j]*E_plr[3][i][j+1] - am[i][j]*E_plr[2][i][j] - (f4W-f4E) ) / dist;

                double drho = UTL::minmod(
                    rho_star - rho_plr[4][i][j],
                    rho_plr[5][i][j+1] - rho_star,
                    rho_star-rho_plr[6][i][j],
                    rho_plr[7][i][j+1]-rho_star
                );

                double dmx = UTL::minmod(
                    mx_star - mx_plr[4][i][j],
                    mx_plr[5][i][j+1] - mx_star,
                    mx_star-mx_plr[6][i][j],
                    mx_plr[7][i][j+1]-mx_star
                );

                double dmy = UTL::minmod(
                    my_star - my_plr[4][i][j],
                    my_plr[5][i][j+1] - my_star,
                    my_star-my_plr[6][i][j],
                    my_plr[7][i][j+1]-my_star
                );

                double dE = UTL::minmod(
                    E_star - E_plr[4][i][j],
                    E_plr[5][i][j+1] - E_star,
                    E_star-E_plr[6][i][j],
                    E_plr[7][i][j+1]-E_star
                );

                F1[i][j] = ( (ap[i][j]*f1E - am[i][j]*f1W) + prod*( rho_plr[3][i][j+1] - rho_plr[2][i][j] - drho ) ) / dist;
                F2[i][j] = ( (ap[i][j]*f2E - am[i][j]*f2W) + prod*( mx_plr[3][i][j+1] - mx_plr[2][i][j] - dmx ) ) / dist;
                F3[i][j] = ( (ap[i][j]*f3E - am[i][j]*f3W) + prod*( my_plr[3][i][j+1] - my_plr[2][i][j] - dmy ) ) / dist;
                F4[i][j] = ( (ap[i][j]*f4E - am[i][j]*f4W) + prod*( E_plr[3][i][j+1] - E_plr[2][i][j] - dE ) ) / dist;
            }
            else{
                F1[i][j] = 0.5*( f1E + f1W );
                F2[i][j] = 0.5*( f2E + f2W );
                F3[i][j] = 0.5*( f3E + f3W );
                F4[i][j] = 0.5*( f4E + f4W );
            }

        }
    }

    //! lspy calculation
    vvd bp(row, vd(col,0));
    vvd bm(row, vd(col,0));

    // CURH Numerical flux
    vvd G1(row, vd(col,0));
    vvd G2(row, vd(col,0));
    vvd G3(row, vd(col,0));
    vvd G4(row, vd(col,0));

    for(int i=0 ; i<row-1 ; i++){
        for(int j=1 ; j<col-1 ; j++){
            // first we need to calculate the required primitive variables
            double uN = mx_plr[0][i][j] / rho_plr[0][i][j];
            double uS = mx_plr[1][i+1][j] / rho_plr[1][i+1][j];
            double vN = my_plr[0][i][j] / rho_plr[0][i][j];
            double vS = my_plr[1][i+1][j] / rho_plr[1][i+1][j];
            double pN = ( CTS::GAMMA - 1 )*( E_plr[0][i][j] - 0.5*rho_plr[0][i][j]*( pow(uN, 2) + pow(vN,2) ) );
            double pS = ( CTS::GAMMA - 1 )*( E_plr[1][i+1][j] - 0.5*rho_plr[1][i+1][j]*( pow(uS, 2) + pow(vS,2) ) );
        
            // now we can calulcate the Local Speeds of propagation
            bp[i][j] = max(
                0.0,
                max(
                    vS + sqrt( CTS::GAMMA*pS/rho_plr[1][i+1][j] ),
                    vN + sqrt( CTS::GAMMA*pN/rho_plr[0][i][j] )
                )
            );

            bm[i][j] = min(
                0.0,
                min(
                    vS - sqrt( CTS::GAMMA*pS/rho_plr[1][i+1][j] ),
                    vN - sqrt( CTS::GAMMA*pN/rho_plr[0][i][j] )
                )
            );

            //! The following optimization on the local speeds of propagation is part of CURH Scheme
            //! if it is excluded then the code runs for CU scheme
            if( CTS::ScaleSw ){
                double differ1 = abs(0.5*( rho_plr[1][i+1][j]*pow(vS,2) - rho_plr[0][i][j]*pow(vN,2) ) + ( (pS-pN) / (CTS::GAMMA-1) ));
                double differ2 = 0.5*abs( rho_plr[1][i+1][j]*pow(uS,2) - rho_plr[0][i][j]*pow(uN,2) );
                double differ = sqrt(pow(differ1,2)+pow(differ2,2));

                double alpha;

                if( differ > CTS::EPSILON ){
                    alpha = differ1/differ;
                }
                else{
                    alpha = 0;
                }

                bp[i][j] = max(
                    0.0,
                    max(
                        vS + alpha*sqrt(CTS::GAMMA*pS/rho_plr[1][i+1][j]),
                        vN + alpha*sqrt(CTS::GAMMA*pN/rho_plr[0][i][j])
                    )
                );

                bm[i][j] = min(
                    0.0,
                    min(
                        vS - alpha*sqrt(CTS::GAMMA*pS/rho_plr[1][i+1][j]),
                        vN - alpha*sqrt(CTS::GAMMA*pN/rho_plr[0][i][j])
                    )
                );
            }

            double g1N = my_plr[0][i][j];
            double g2N = my_plr[0][i][j]*uN;
            double g3N = my_plr[0][i][j]*vN + pN;
            double g4N = vN*( E_plr[0][i][j] + pN );

            double g1S = my_plr[1][i+1][j];
            double g2S = my_plr[1][i+1][j]*uS;
            double g3S = my_plr[1][i+1][j]*vS + pS;
            double g4S = vS*( E_plr[1][i+1][j] + pS );

            if( CTS::ShrinkSw ){
                double dg1 = g1S - g1N;
                double dg2 = g2S - g2N;
                double dg3 = g3S - g3N;
                double dg4 = g4S - g4N;

                double du1 = rho_plr[1][i+1][j] - rho_plr[0][i][j];
                double du2 = mx_plr[1][i+1][j] - mx_plr[0][i][j];
                double du3 = my_plr[1][i+1][j] - my_plr[0][i][j];
                double du4 = E_plr[1][i+1][j] - E_plr[0][i][j];

                double du1eps, du2eps, du3eps, du4eps;

                if( du1 > 0 ){
                    du1eps = max(du1, 1.0E-10);
                }
                else{  
                    du1eps = min(du1, -1.0E-10);
                }

                if( du2 > 0 ){
                    du2eps = max(du2, 1.0E-10);
                }
                else{  
                    du2eps = min(du2, -1.0E-10);
                }

                if( du3 > 0 ){
                    du3eps = max(du3, 1.0E-10);
                }
                else{  
                    du3eps = min(du3, -1.0E-10);
                }

                if( du4 > 0 ){
                    du4eps = max(du4, 1.0E-10);
                }
                else{  
                    du4eps = min(du4, -1.0E-10);
                }

                double ratio1 = 2.0*dg1 / (du1+du1eps);
                double ratio2 = 2.0*dg2 / (du2+du2eps);
                double ratio3 = 2.0*dg3 / (du3+du3eps);
                double ratio4 = 2.0*dg4 / (du4+du4eps);

                double ratiomin = min(ratio1, min(ratio2, min(ratio3, ratio4)));
                double ratiomax = max(ratio1, max(ratio2, max(ratio3, ratio4)));

                if( ratiomax > 0 ){
                    bp[i][j] = min( bp[i][j] , ratiomax );
                    bm[i][j] = max( bm[i][j] , -ratiomax );
                }

                if( ratiomin < 0 ){
                    bp[i][j] = min( bp[i][j] , -ratiomin );
                    bm[i][j] = max( bm[i][j] , ratiomin );
                }
            }

            // CURH Numerical flux calculation starts
            double dist = bp[i][j] - bm[i][j];
            double prod = bp[i][j]*bm[i][j];

            if( dist > CTS::EPSILON ){
                double rho_star = ( bp[i][j]*rho_plr[1][i+1][j] - bm[i][j]*rho_plr[0][i][j] - (g1S-g1N) ) / dist;
                double mx_star = ( bp[i][j]*mx_plr[1][i+1][j] - bm[i][j]*mx_plr[0][i][j] - (g2S-g2N) ) / dist;
                double my_star = ( bp[i][j]*my_plr[1][i+1][j] - bm[i][j]*my_plr[0][i][j] - (g3S-g3N) ) / dist;
                double E_star = ( bp[i][j]*E_plr[1][i+1][j] - bm[i][j]*E_plr[0][i][j] - (g4S-g4N) ) / dist;

                double drho = UTL::minmod(
                    rho_star - rho_plr[4][i][j],
                    rho_plr[6][i+1][j] - rho_star,
                    rho_star - rho_plr[5][i][j],
                    rho_plr[7][i+1][j] - rho_star
                );

                double dmx = UTL::minmod(
                    mx_star - mx_plr[4][i][j],
                    mx_plr[6][i+1][j] - mx_star,
                    mx_star - mx_plr[5][i][j],
                    mx_plr[7][i+1][j] - mx_star
                );

                double dmy = UTL::minmod(
                    my_star - my_plr[4][i][j],
                    my_plr[6][i+1][j] - my_star,
                    my_star - my_plr[5][i][j],
                    my_plr[7][i+1][j] - my_star
                );

                double dE = UTL::minmod(
                    E_star - E_plr[4][i][j],
                    E_plr[6][i+1][j] - E_star,
                    E_star - E_plr[5][i][j],
                    E_plr[7][i+1][j] - E_star
                );

                G1[i][j] = ( (bp[i][j]*g1N - bm[i][j]*g1S) + prod*( ( rho_plr[1][i+1][j] - rho_plr[0][i][j] ) - drho ) ) / dist;
                G2[i][j] = ( (bp[i][j]*g2N - bm[i][j]*g2S) + prod*( ( mx_plr[1][i+1][j] - mx_plr[0][i][j] ) - dmx ) ) / dist;
                G3[i][j] = ( (bp[i][j]*g3N - bm[i][j]*g3S) + prod*( ( my_plr[1][i+1][j] - my_plr[0][i][j] ) - dmy ) ) / dist;
                G4[i][j] = ( (bp[i][j]*g4N - bm[i][j]*g4S) + prod*( ( E_plr[1][i+1][j] - E_plr[0][i][j] ) - dE ) ) / dist;
            }  
            else{
                G1[i][j] = 0.5*( g1N + g1S );                
                G2[i][j] = 0.5*( g2N + g2S );                
                G3[i][j] = 0.5*( g3N + g3S );                
                G4[i][j] = 0.5*( g4N + g4S );                
            }

        }
    }

    vvvd F = {F1, F2, F3, F4};
    vvvd G = {G1, G2, G3, G4};

    if( update_time ){
        double amax=0.0;
        double bmax=0.0;

        for(int i=1 ; i<row-1 ; i++){
            for(int j=0 ; j<col-1 ; j++){
                amax = max(
                    amax,
                    max(
                        ap[i][j],
                        -am[i][j]
                    )
                );
            }
        }

        for(int i=0 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                bmax = max(
                    bmax,
                    max(
                        bp[i][j],
                        -bm[i][j]
                    )
                );
            }
        }

        dt = CTS::CFL * min(dx/amax , dy/bmax);

        if( t+dt > time.second ){
            dt = time.second - t;
        }
    }  

    return {F,G};
}

// void CURH::update_dt(vvvd& lspx, vvvd& lspy){

//     int row = lspx[0].size();
//     int col = lspx[0][0].size();

//     double amax=0.0;
//     double bmax=0.0;

//     for(int i=1 ; i<row-1 ; i++){
//         for(int j=0 ; j<col-1 ; j++){
//             amax = max(
//                 amax,
//                 max(
//                     lspx[0][i][j],
//                     -lspx[1][i][j]
//                 )
//             );
//         }
//     }

//     for(int i=0 ; i<row-1 ; i++){
//         for(int j=1 ; j<col-1 ; j++){
//             bmax = max(
//                 bmax,
//                 max(
//                     lspy[0][i][j],
//                     -lspy[1][i][j]
//                 )
//             );
//         }
//     }

//     dt = CTS::CFL * min(dx/amax , dy/bmax);

//     if( t+dt > time.second ){
//         dt = time.second - t;
//     }

// }

void CURH::update_conserved_variables(){

    // get the curh flux corresponding to conserved variables at time t=t0
    vvvvd curh_flux_1 = curh_flux(rho, mx, my, E, true);

    vvvd F = curh_flux_1[0];
    vvvd G = curh_flux_1[1];

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
        for(int j=1 ; j<col-1 ; j++){
            rho1[i-1][j-1] = rho[i][j] - ( LAMBDA*(F[0][i][j] - F[0][i][j-1]) + MU*(G[0][i][j] - G[0][i-1][j]) );
            mx1[i-1][j-1] = mx[i][j] - ( LAMBDA*(F[1][i][j] - F[1][i][j-1]) + MU*(G[1][i][j] - G[1][i-1][j]) );
            my1[i-1][j-1] = my[i][j] - ( LAMBDA*(F[2][i][j] - F[2][i][j-1]) + MU*(G[2][i][j] - G[2][i-1][j]) );
            E1[i-1][j-1] = E[i][j] - ( LAMBDA*(F[3][i][j] - F[3][i][j-1]) + MU*(G[3][i][j] - G[3][i-1][j]) );
        }
    }

    EXC::extend_matrix(rho1, bc, "Density");
    EXC::extend_matrix(mx1, bc, "MomentumX");
    EXC::extend_matrix(my1, bc, "MomentumY");
    EXC::extend_matrix(E1, bc, "Energy");

    vvvvd curh_flux_2 = curh_flux(rho1, mx1, my1, E1, false);

    F = curh_flux_2[0];
    G = curh_flux_2[1];

    // STAGE-2
    vvd rho2(row-2, vd(col-2));
    vvd mx2(row-2, vd(col-2));
    vvd my2(row-2, vd(col-2));
    vvd E2(row-2, vd(col-2));

    for(int i=1 ; i<row-1 ; i++){
        for(int j=1 ; j<col-1 ; j++){
            rho2[i-1][j-1] = ( (3.0*rho[i][j] + rho1[i][j]) - ( LAMBDA*(F[0][i][j] - F[0][i][j-1]) + MU*(G[0][i][j] - G[0][i-1][j]) ) ) / 4.0;
            mx2[i-1][j-1] = ( (3.0*mx[i][j] + mx1[i][j]) - ( LAMBDA*(F[1][i][j] - F[1][i][j-1]) + MU*(G[1][i][j] - G[1][i-1][j]) ) ) / 4.0;
            my2[i-1][j-1] = ( (3.0*my[i][j] + my1[i][j]) - ( LAMBDA*(F[2][i][j] - F[2][i][j-1]) + MU*(G[2][i][j] - G[2][i-1][j]) ) ) / 4.0;
            E2[i-1][j-1] = ( (3.0*E[i][j] + E1[i][j]) - ( LAMBDA*(F[3][i][j] - F[3][i][j-1]) + MU*(G[3][i][j] - G[3][i-1][j]) ) ) / 4.0;
        }
    }

    EXC::extend_matrix(rho2, bc, "Density");
    EXC::extend_matrix(mx2, bc, "MomentumX");
    EXC::extend_matrix(my2, bc, "MomentumY");
    EXC::extend_matrix(E2, bc, "Energy");

    vvvvd curh_flux_3 = curh_flux(rho2, mx2, my2, E2, false);

    F = curh_flux_3[0];
    G = curh_flux_3[1];

    // STAGE-3
    vvd rhon(row-2, vd(col-2));
    vvd mxn(row-2, vd(col-2));
    vvd myn(row-2, vd(col-2));
    vvd En(row-2, vd(col-2));

    for(int i=1 ; i<row-1 ; i++){
        for(int j=1 ; j<col-1 ; j++){
            rhon[i-1][j-1] = ( ( rho[i][j] + 2.0*rho2[i][j] ) - 2.0*( LAMBDA*(F[0][i][j] - F[0][i][j-1]) + MU*(G[0][i][j] - G[0][i-1][j]) ) ) / 3.0; 
            mxn[i-1][j-1] = ( ( mx[i][j] + 2.0*mx2[i][j] ) - 2.0*( LAMBDA*(F[1][i][j] - F[1][i][j-1]) + MU*(G[1][i][j] - G[1][i-1][j]) ) ) / 3.0; 
            myn[i-1][j-1] = ( ( my[i][j] + 2.0*my2[i][j] ) - 2.0*( LAMBDA*(F[2][i][j] - F[2][i][j-1]) + MU*(G[2][i][j] - G[2][i-1][j]) ) ) / 3.0; 
            En[i-1][j-1] = ( ( E[i][j] + 2.0*E2[i][j] ) - 2.0*( LAMBDA*(F[3][i][j] - F[3][i][j-1]) + MU*(G[3][i][j] - G[3][i-1][j]) ) ) / 3.0; 
        }
    }

    EXC::extend_matrix(rhon, bc, "Density");
    EXC::extend_matrix(mxn, bc, "MomentumX");
    EXC::extend_matrix(myn, bc, "MomentumY");
    EXC::extend_matrix(En, bc, "Energy");

    //! update the conserved variables
    rho = rhon;
    mx = mxn;
    my = myn;
    E = En;

}


