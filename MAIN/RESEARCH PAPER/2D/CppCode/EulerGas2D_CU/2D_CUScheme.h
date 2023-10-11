/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<unordered_map>
#include<cmath>

#include "Constants.h"
#include "Utility.h"
#include "ExtendCells.h"
#include "Flux.h"
#include "PrimitiveVariables.h"
#include "OutputResult.h"

using namespace std;

namespace CTS = Constants;
namespace CTR = ConstantsRT;
namespace UTL = Utility;
namespace OPR = OutputResult;
namespace EXC = ExtendCells;
namespace PRV = PrimitiveVariables;

#define ll long long int
#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define vvvvd vector< vector< vector< vector<double> > > >
#define pdd pair<double,double>
#define pii pair<ll,ll>
#define umap_ss unordered_map<string,string>

class EulerSystem_2D{

// class variables represent the variables of the Riemann problem
// All the variables are initialized in the class constructor
private:

    string mode; // CU ; REF
    string problem; // MCW ; 2DR ; EXP ; IMP ; KHI ; RTI

    pdd domainX; // computational domain in the x-direction
    pdd domainY; // computational domain in the y-direction

    pii Nxy; // number of cells in the x and y directions

    pdd time; // initial and final time

    double t, dt; // current time and time step
    double dx, dy; // length of individual cell in the Finite Volume Grid in x and y directions

    umap_ss bc; // boundary conditions on the 4 boundaries [N, S, E, W]

    vd gridX; // grid points in x-direction
    vd gridY; // grid points in y-direction

    // 2D vectors for the Conserved variables which are evolved in time using the CU Scheme
    vvd rho; // density
    vvd mx; // x-momentum
    vvd my; // y-momentum
    vvd E; // Energy

public:

    // Function to run the complete CU scheme and store data at every 5th iteration
    // used for animating the results
    void RunCU_complete(){
        
    }

    // Function to run the partial CU scheme and store data only at start and end of the iterations
    // used for static image plots initially and finally
    void RunCU_partial(){

    }

    /*
        ! STEP-1 Initialize the Initial conditions

        Functions to initialize the variables and create and export the environment
        are defined as follows
    */

    /*
        Constructor takes the problem string as input to initialize the variables
    
        -> [MCW] Moving Contact Wave
        -> [2DR] 2-Dimensional Riemann Problem
        -> [EXP] Explosion Problem
        -> [IMP] Implosion Problem
        -> [KHI] Kelvin Helmholtz Instability
        -> [RTI] Raleigh Taylor Instability
    */
    EulerSystem_2D(string mode, string problem) : mode(mode), problem(problem){
        // log message
        cout << "---LOG--- Running CU Scheme for - " << problem << "[" << mode << "] --- " << endl;

        // initialize the environment and set up initial variables and matrices
        
        initialize_environment(); // initialize environment variables
        make_grids(); // make one dimensional grids
        initialize_riemann_problem(); // initialize the conserved variables

        // write the one dimensional grids to the external file
        OPR::write_grids(mode, gridX, gridY);

        // export environment variables
        export_env();
    }

    // export environment variables for use in plotting
    void export_env(){
        cout << "---LOG--- Exporting environment variables --- " << endl;

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
        ss << domainX.first << " " << domainX.second << endl;
        ss << domainY.first << " " << domainY.second << endl;
        ss << time.first << " " << time.second << endl;
        ss << Nxy.first << endl;
        ss << Nxy.second << endl;

        fout << ss.str() << endl;

        fout.close();
    }

    // Function to initialize the variables according to the Problem being used
    void initialize_environment(){
        cout << "---LOG--- Initializing problem environment --- " << endl;
        
        if( mode == "CU" ){ // CU PLot with less precise initial data
        
            if( problem == "MCW" ){ // Moving Contact Wave
                domainX.first = -0.2;
                domainX.second = 0.2;

                domainY.first = 0.0;
                domainY.second = 0.8;

                dx = (double)1/250;
                dy = (double)1/250;

                // calculate Nx and Ny for
                Nxy.first = (ll)( (domainX.second-domainX.first) / dx );
                Nxy.second = (ll)( (domainY.second-domainY.first) / dy );

                time.first = 0.0;
                time.second = 2.0;

                t = time.first;

                // free conditions
                bc["N"] = "FREE";
                bc["S"] = "FREE";
                bc["E"] = "FREE";
                bc["W"] = "FREE";
            }
            else if( problem == "2DR" ){ // 2-Dimensional Riemann Problem

            }
            else if( problem == "EXP" ){ // Explosion Problem

            }
            else if( problem == "IMP" ){ // Implosion Problem

            }
            else if( problem == "KHI" ){ // Kelvin Helmholtz Instability

            }
            else if( problem == "RTI" ){ // Raleigh Taylor Instability

            }
            else{
                cout << "---ERROR--- Please enter correct riemann problem --- " << endl;
            }

        }
        //! computational grid is 4x more fine in Reference plots
        else if( mode == "REF" ){ // REF Plot with more precise initial data

            if( problem == "MCW" ){ // Moving Contact Wave
                domainX.first = -0.2;
                domainX.second = 0.2;

                domainY.first = 0.0;
                domainY.second = 0.8;

                dx = (double)1/1000; 
                dy = (double)1/1000;

                // calculate Nx and Ny for
                Nxy.first = (ll)( (domainX.second-domainX.first) / dx );
                Nxy.second = (ll)( (domainY.second-domainY.first) / dy );

                time.first = 0.0;
                time.second = 2.0;

                t = time.first;

                // free conditions
                bc["N"] = "FREE";
                bc["S"] = "FREE";
                bc["E"] = "FREE";
                bc["W"] = "FREE";
            }
            else if( problem == "2DR" ){ // 2-Dimensional Riemann Problem

            }
            else if( problem == "EXP" ){ // Explosion Problem

            }
            else if( problem == "IMP" ){ // Implosion Problem

            }
            else if( problem == "KHI" ){ // Kelvin Helmholtz Instability

            }
            else if( problem == "RTI" ){ // Raleigh Taylor Instability

            }
            else{
                cout << "---ERROR--- Please enter correct riemann problem --- " << endl;
            }

        }
        else{
            cout << "---ERROR--- Please enter a correct mode ---" << endl;
        }

    }

    // Function to make the one dimensional grids in x-direction and y-direction
    void make_grids(){
        cout << "---LOG--- Initializing one dimensional grids ---" << endl;

        for(int i=0 ; i<Nxy.first ; i++){
            gridX.push_back( domainX.first + dx*(i+0.5) );
        }
        
        for(int i=0 ; i<Nxy.second ; i++){
            gridY.push_back( domainY.first + dy*(i+0.5) );
        }

    }

    // Function to Initialize the Conserved Vectors which are class variables
    // based on the given problems
    void initialize_riemann_problem(){
        cout << "---LOG--- Initializing conserved variables ---" << endl;

        if( problem == "MCW" ){ // Moving Contact Wave

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
        else if( problem == "2DR" ){ // 2-Dimensional Riemann Problem

        }
        else if( problem == "EXP" ){ // Explosion Problem

        }
        else if( problem == "IMP" ){ // Implosion Problem

        }
        else if( problem == "KHI" ){ // Kelvin Helmholtz Instability

        }
        else if( problem == "RTI" ){ // Raleigh Taylor Instability

        }
        else{
            cout << "---ERROR--- Please enter correct problem ---" << endl;
        }

        // Extend cells for the initial problem
        EXC::extend_matrix(rho, bc, "Density");
        EXC::extend_matrix(rho, bc, "MomentumX");
        EXC::extend_matrix(rho, bc, "MomentumY");
        EXC::extend_matrix(rho, bc, "Energy");

    }

    /*
        ! STEP-2 Run the CU scheme
        Main functions to run the CU scheme are defined as follows
    */

    //! -> Construct the Piecewise Linear Reconstruction

    // Function to calculate the plr
    // [N, S, E, W, NE, NW, SE, SW]
    vvvd get_plr(vvd& var, string var_name){
        // first get the slopes
        vvd slx = get_slx(var);
        vvd sly = get_sly(var);
    
        int row = var.size();
        int col = var[0].size();

        vvvd plr;

        vvd n(row-2, vd(col-2));
        vvd s(row-2, vd(col-2));
        vvd e(row-2, vd(col-2));
        vvd w(row-2, vd(col-2));
        vvd ne(row-2, vd(col-2));
        vvd nw(row-2, vd(col-2));
        vvd se(row-2, vd(col-2));
        vvd sw(row-2, vd(col-2));
    
        for(int i=1 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                n[i-1][j-1] = var[i][j] + 0.5*sly[i][j];
                s[i-1][j-1] = var[i][j] - 0.5*sly[i][j];
                e[i-1][j-1] = var[i][j] + 0.5*slx[i][j];
                w[i-1][j-1] = var[i][j] - 0.5*slx[i][j];
                ne[i-1][j-1] = var[i][j] + 0.5*slx[i][j] + 0.5*sly[i][j];
                nw[i-1][j-1] = var[i][j] - 0.5*slx[i][j] + 0.5*sly[i][j];
                se[i-1][j-1] = var[i][j] + 0.5*slx[i][j] - 0.5*sly[i][j];
                sw[i-1][j-1] = var[i][j] - 0.5*slx[i][j] - 0.5*sly[i][j]; 
            }
        }

        EXC::extend_matrix(n, bc, var_name);
        EXC::extend_matrix(s, bc, var_name);
        EXC::extend_matrix(e, bc, var_name);
        EXC::extend_matrix(w, bc, var_name);
        EXC::extend_matrix(ne, bc, var_name);
        EXC::extend_matrix(nw, bc, var_name);
        EXC::extend_matrix(se, bc, var_name);
        EXC::extend_matrix(sw, bc, var_name);

        plr.push_back(n);
        plr.push_back(s);
        plr.push_back(e);
        plr.push_back(w);
        plr.push_back(ne);
        plr.push_back(nw);
        plr.push_back(se);
        plr.push_back(sw);

        return plr;

    }

    // Function to calculate the slx
    vvd get_slx(vvd& var){
        int row = var.size();
        int col = var[0].size();

        vvd slx(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                double v1 = CTS::THETA*( var[i][j] - var[i][j-1] );
                double v2 = 0.5*( var[i][j+1] - var[i][j-1] );
                double v3 = CTS::THETA*( var[i][j+1] - var[i][j] );

                slx[i-1][j-1] = UTL::minmod( v1 , v2 , v3 );
            }
        }

        EXC::extend_matrix(slx, bc, "SLOPES");

        return slx;
    }        

    // Function to calculate the sly
    vvd get_sly(vvd& var){
        int row = var.size();
        int col = var[0].size();

        vvd sly(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                double v1 = CTS::THETA*( var[i][j] - var[i-1][j] );
                double v2 = 0.5*( var[i+1][j] - var[i-1][j] );
                double v3 = CTS::THETA*( var[i+1][j] - var[i][j] );

                sly[i-1][j-1] = UTL::minmod( v1 , v2 , v3 );
            }
        }

        EXC::extend_matrix(sly, bc, "SLOPES");

        return sly;
    }  

    //! -> Calculate the One Sided Local Speeds of Propagation

    // Function to compute the One sided local speeds of propagation in the x-direction
    vvvd get_lsp_x(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr){
        double GAMMA;
        
        if( problem == "RTI" ){
            GAMMA = CTR::GAMMA;
        }
        else{
            GAMMA = CTS::GAMMA;
        }

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

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){
                ap[i][j] = max(
                    0.0 ,
                    max(
                        uW[i][j+1] + sqrtf( abs( (GAMMA*pW[i][j+1])/rho_plr[3][i][j+1] ) ),
                        uE[i][j] + sqrtf( abs( (GAMMA*pE[i][j]) / rho_plr[2][i][j] ) )
                    )
                );

                am[i][j] = min(
                    0.0,
                    min(
                        uW[i][j+1] - sqrtf( abs( (GAMMA*pW[i][j+1])/rho_plr[3][i][j+1] ) ),
                        uE[i][j] - sqrtf( abs( (GAMMA*pE[i][j])/rho_plr[2][i][j] ) )
                    )
                );
            }
        }

        lsp.push_back(ap);
        lsp.push_back(am);

        return lsp;
    }

    // Function to compute the One sided local speeds of propagation in the y-direction
    vvvd get_lsp_y(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr){
        double GAMMA;
        
        if( problem == "RTI" ){
            GAMMA = CTR::GAMMA;
        }
        else{
            GAMMA = CTS::GAMMA;
        }

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

        vvd bp(row, vd(col));
        vvd bm(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){
                bp[i][j] = max(
                    0.0 ,
                    max(
                        vW[i+1][j] + sqrtf( abs((GAMMA*pW[i+1][j])/rho_plr[3][i+1][j]) ),
                        vE[i][j] + sqrtf( abs((GAMMA*pE[i][j]) / rho_plr[2][i][j]) )
                    )
                );

                bm[i][j] = min(
                    0.0,
                    min(
                        vW[i+1][j] - sqrtf( abs((GAMMA*pW[i+1][j])/rho_plr[3][i+1][j]) ),
                        vE[i][j] - sqrtf( abs((GAMMA*pE[i][j])/rho_plr[2][i][j]) )
                    )
                );
            }
        }

        lsp.push_back(bp);
        lsp.push_back(bm);

        return lsp;
    }

    //! -> Calculate the Built in Anti-Diffusion term

    //! -> Calculate the CU Numerical Flux

    //! -> Update the values of Conserved variables using 3-stage SSPRK Scheme

};
