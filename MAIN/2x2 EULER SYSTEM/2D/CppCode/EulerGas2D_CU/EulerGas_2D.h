/*
    * @author - lakshya Bamne (20MA20029)
    * @supervisor - Prof. Naveen Kumar Garg

    ! Implementation of CU scheme for modified Flux functions in Euler's Equations
*/

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<unordered_map>

#include "ExtendCells.h"
#include "Constants.h"
#include "Utility.h"
#include "Flux.h"
#include "OutputResult.h"

using namespace std;

#define ll long long int
#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define vvvvd vector< vector< vector< vector<double> > > >
#define umap_ss unordered_map<string,string>
#define pdd pair<double,double>
#define pii pair<int,int>
#define pss pair<string,string>

namespace EXC = ExtendCells;
namespace CTS = Constants;
namespace UTL = Utility;
namespace FLX = Flux;
namespace OPR = OutputResult;

class ModifiedEulerSystem_CU{
private:
    int shock_type;
    
    // initial velocity values for the riemann problem
    vd state;
    
    // density value for the problem
    double rho_val;

    // domain end points in x and y directions
    pdd domainX;
    pdd domainY;

    // time start and end values
    pdd time;

    // number of grid cells in the domain
    pii grid_cells;

    // boundary conditions are neumann(free) for all boundaries
    umap_ss bc;

    // variables required to make the Initial Grid for the problem
    double dx, dy;
    double t, dt;

    // Riemann Grids
    vd gridX, gridY;

    vvd u;
    vvd rho;

public:
    /*
        ! Functions to run the program in different modes (Complete, Partial)
    */
    void RunCU_complete(){

        while( t < time.second ){
            cout << "t = " << t << " | dt = " << dt << endl;
            
            //! write variable values at each step for making animations
            OPR::write_matrix(this->shock_type, "VELOCITY", u);
            OPR::write_matrix(this->shock_type, "DENSITY", rho);

            // function to update the conservative variables internally
            // also update dt for next iteration
            get_next_vars();

            // update new time step calculated using CFL conditions
            t += dt;

        }
    }

    void RunCU_partial(){
        //! write initial values for plotting
        OPR::write_matrix(this->shock_type, "VELOCITY", u);
        OPR::write_matrix(this->shock_type, "DENSITY", rho);

        while( t < time.second ){
            cout << "------------------------------------------------------------------------------------" << endl;
            cout << "t = " << t << " | dt = " << dt << endl;
            
            // function to update the conservative variables internally
            // also update dt for next iteration
            get_next_vars();

            // update new time step calculated using CFL conditions
            t += dt;

        }

        //! write final values for plotting
        OPR::write_matrix(this->shock_type, "VELOCITY", u);
        OPR::write_matrix(this->shock_type, "DENSITY", rho);
    }

    // function to export the environment variables for the problem to be used during plotting
    void export_env(int shock_type){
        stringstream ss1;
        ss1 << "env/Environment" << shock_type << ".txt"; 
        
        string file_name = ss1.str();
        
        stringstream ss;
        
        ofstream fout;
        fout.open(file_name , ios::app);

        string mode;

        switch(shock_type){
            case 1: // JS-JS-SJ
                mode = "JS-JS-SJ";
                
                break;
            case 2: // JR-JS-SJ
                mode = "JR-JS-SJ";
                
                break;
            case 3: // RJ-SJ-JS
                mode = "RJ-SJ-JS";
                
                break;
            case 4: // JR-JS-RJ
                mode = "JR-JS-RJ";
                
                break;
            case 5: // R-JS-R
                mode = "R-JS-R";
                
                break;
            case 6: // RJ-SJ-JR
                mode = "RJ-SJ-JR";
                
                break;
            case 7: // JR-JR-RJ
                mode = "JR-JR-RJ";
                
                break;
            case 8: // R-JR-R
                mode = "R-JR-R";
                
                break;
            case 9: // SD-R-RJ
                mode = "SD-R-RJ";
                
                break;
            case 10: // SD-RJ-SJ
                mode = "SD-RJ-SJ";
                
                break;
            case 11: //  SD-RJ-SD
                mode = "SD-RJ-SD";
                
                break;
            case 12: // SD-SJ-SD
                mode = "SD-SJ-SD";
                break;
            
            default:
                cout << "---ERROR--- Please enter a valid initial state ---" << endl;
                break;
        }

        ss << mode << endl;
        ss << domainX.first << " " << domainX.second << endl;
        ss << domainY.first << " " << domainY.second << endl;
        ss << time.first << " " << time.second << endl;
        ss << grid_cells.first << endl;
        ss << grid_cells.second << endl;

        fout << ss.str() << endl;

        fout.close();
    }

    /*
        ! STEP-1 : Initialize the variables and the computational grid
    */

    // Class constructor is used to initialize the Riemann Problem
    ModifiedEulerSystem_CU(int shock_type){

        // initializing the type of shock in the class variable for later usage
        this->shock_type = shock_type;

        rho_val = 0.77;
        domainX = {-1,1};
        domainY = {-1,1};
        time = {0,0.25};
        grid_cells = {600,600};
        // grid_cells = {10,10};

        bc["N"] = "FREE";
        bc["S"] = "FREE";
        bc["E"] = "FREE";
        bc["W"] = "FREE";

        dx = (domainX.second - domainX.first)/grid_cells.first;
        dy = (domainY.second - domainY.first)/grid_cells.second;

        t = time.first;

        // initializing the variables according to the problem
        switch(shock_type){
            case 1: // JS-JS-SJ
                state.push_back(-0.56);            
                state.push_back(-0.37);            
                state.push_back(-0.15);

                break;
            case 2: // JR-JS-SJ
                state.push_back(-0.37);            
                state.push_back(-0.56);            
                state.push_back(-0.15);

                break;
            case 3: // RJ-SJ-JS
                state.push_back(0.37);            
                state.push_back(0.15);
                state.push_back(0.56);            

                break;
            case 4: // JR-JS-RJ
                state.push_back(-0.17);
                state.push_back(-0.56);            
                state.push_back(-0.37);            

                break;
            case 5: // R-JS-R
                state.push_back(0.17);
                state.push_back(-0.56);            
                state.push_back(-0.37);            

                break;
            case 6: // RJ-SJ-JR
                state.push_back(0.56);            
                state.push_back(0.15);
                state.push_back(0.37);            

                break;
            case 7: // JR-JR-RJ
                state.push_back(-0.12);            
                state.push_back(-0.37);            
                state.push_back(-0.56);

                break;
            case 8: // R-JR-R
                state.push_back(0.15);            
                state.push_back(-0.37);            
                state.push_back(-0.56);

                break;
            case 9: // SD-R-RJ
                state.push_back(-0.15);            
                state.push_back(0.45);            
                state.push_back(-0.37);

                break;
            case 10: // SD-RJ-SJ
                state.push_back(-0.37);            
                state.push_back(0.25);            
                state.push_back(-0.15);

                break;
            case 11: //  SD-RJ-SD
                state.push_back(-0.37);            
                state.push_back(0.37);            
                state.push_back(0.15);
                
                break;
            case 12: // SD-SJ-SD
                state.push_back(-0.56);            
                state.push_back(0.15);            
                state.push_back(0.37);

                break;
            
            default:
                cout << "---ERROR--- Please enter a valid initial state ---" << endl;
                break;
        }

        cout << "---LOG--- Successfully initialized variables for the Riemann Problem ---" << endl;

        // after initializing the variables, we have to initialize the grid at t = 0
        InitializeGrid();

        // export the environment variables for use in plotting
        export_env(this->shock_type);

        // write the initial grid to output file for plotting
        OPR::write_grids(this->shock_type, gridX, gridY);

        cout << "---LOG--- Initialized Riemann Grid ---" << endl;
    }

    // function to show the riemann grid created
    void show_cons_vars(){
        cout << "----------------------------------VARIABLES---------------------------------------" << endl;

        cout << "Velocity(u) :-" << endl;
        for(int j=u.size()-1 ; j>=0 ; j--){
            for(auto i : u[j]){
                cout << i << " ";
            }
            cout << endl;
        }

        cout << "Density(rho) :-" << endl;
        for(int j=rho.size()-1 ; j>=0 ; j--){
            for(auto i : rho[j]){
                cout << i << " ";
            }
            cout << endl;
        }

        cout << "----------------------------------------------------------------------------------" << endl;
    }

    // function to initialize the Riemann Grid for the problem
    void InitializeGrid(){
        // first initialize the x and y vectors using which 2D grid is created
        for(int i=0 ; i<grid_cells.first ; i++){
            gridX.push_back( domainX.first + (i+0.5)*( dx ) );
        }

        for(int i=0 ; i<grid_cells.second ; i++){
            gridY.push_back( domainY.first + (i+0.5)*( dy ) );
        }

        // using the x and y grids we have to initialize the 2D Riemann Grid at t=0
        for(auto j : gridY){
            vd u1;
            vd rho1;

            for(auto i : gridX){
                rho1.push_back(rho_val);

                if( i>=0 && j>=0 ){
                    u1.push_back(state[0]);
                }
                else if( i<0 && j>0 ){
                    u1.push_back(state[1]);
                }
                else{
                    u1.push_back(state[2]);
                }
            }

            u.push_back(u1);
            rho.push_back(rho1);

        }

        // extend the grid
        ExtendVarMatrices();

    }

    /*
        ! STEP-2 : Run CU Scheme on the variables
    */

    // function to update the variables using the SSPRK scheme and the CU Numerical Fluxes
    void get_next_vars(){
        // using class variables as input
        vvvvd cu_flux = get_cu_flux(u, rho);
        
        vvvd F = cu_flux[0];
        vvvd G = cu_flux[1];

        const double old_dt = dt; // used because dt is changed later which is undesirable

        double LAMBDA = dt/dx;
        double MU = dt/dy;

        int row = rho.size();
        int col = rho[0].size();

        // STAGE-1
        vvd u1(row-2, vd(col-2));
        vvd rho1(row-2, vd(col-2));
    
        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                u1[i-1][j-1] = u[i][j] - ( LAMBDA*( F[0][i][j] - F[0][i][j-1] ) + MU*( G[0][i][j] - G[0][i-1][j] ) );
                rho1[i-1][j-1] = rho[i][j] - ( LAMBDA*( F[1][i][j] - F[1][i][j-1] ) + MU*( G[1][i][j] - G[1][i-1][j] ) );
            }
        }

        EXC::extend_matrix(u1, bc, "VELOCITY");
        EXC::extend_matrix(rho1, bc, "DENSITY");

        // STAGE-2
        vvvvd cu_flux_1 = get_cu_flux(u1, rho1);

        dt = old_dt; // restore value of dt to old value

        vvvd F1 = cu_flux_1[0];
        vvvd G1 = cu_flux_1[1];

        vvd u2(row-2, vd(col-2));
        vvd rho2(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                u2[i-1][j-1] = ( ( 3*u[i][j] + u1[i][j] ) - ( LAMBDA*( F1[0][i][j] - F1[0][i][j-1] ) + MU*( G1[0][i][j] - G1[0][i-1][j] ) ) ) / 4.0;
                rho2[i-1][j-1] = ( ( 3*rho[i][j] + rho1[i][j] ) - ( LAMBDA*( F1[1][i][j] - F1[1][i][j-1] ) + MU*( G1[1][i][j] - G1[1][i-1][j] ) ) ) / 4.0;
            }
        }

        EXC::extend_matrix(u2, bc, "VELOCITY");
        EXC::extend_matrix(rho2, bc, "DENSITY");
        
        // STAGE-3
        vvvvd cu_flux_new = get_cu_flux(u2, rho2);

        dt = old_dt; // restore value of dt to old value

        vvvd FN = cu_flux_new[0];
        vvvd GN = cu_flux_new[1];

        vvd un(row-2, vd(col-2));
        vvd rhon(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                un[i-1][j-1] = ( ( u[i][j] + 2*u2[i][j] ) - 2*( LAMBDA*( FN[0][i][j] - FN[0][i][j-1] ) + MU*( GN[0][i][j] - GN[0][i-1][j] ) ) ) / 3.0;
                rhon[i-1][j-1] = ( ( rho[i][j] + 2*rho2[i][j] ) - 2*( LAMBDA*( FN[1][i][j] - FN[1][i][j-1] ) + MU*( GN[1][i][j] - GN[1][i-1][j] ) ) ) / 3.0;
            }
        }

        EXC::extend_matrix(un, bc, "VELOCITY");
        EXC::extend_matrix(rhon, bc, "DENSITY");

        // Updation of conserved variables to new values
        u = un;
        rho = rhon;
    }

    //! MAIN function to run the CU scheme given a matrix for the variables
    vvvvd get_cu_flux(vvd& u, vvd& rho){
        // get the piecewise linear reconstruction
        vvvd u_plr = get_plr(u, "VELOCITY");
        vvvd rho_plr = get_plr(rho, "DENSITY");

        // get the Local speeds of propagation in x and y directions
        vvvd lsp_x = get_lsp_a(u_plr, rho_plr);
        vvvd lsp_y = get_lsp_b(u_plr, rho_plr);

        // update the time step used in SSPRK
        dt = update_dt(lsp_x, lsp_y);
        
        // we need Flux vectors from the Euler equations later
        vvvd flux_E = FLX::fflux(u_plr[2],rho_plr[2]);
        vvvd flux_W = FLX::fflux(u_plr[3], rho_plr[3]);
        vvvd flux_N = FLX::gflux(u_plr[0],rho_plr[0]);
        vvvd flux_S = FLX::gflux(u_plr[1], rho_plr[1]);

        // now we need to calculate the anti diffusion terms
        vvvd adt_x = get_adt_x(u_plr, rho_plr, lsp_x, flux_E, flux_W);
        vvvd adt_y = get_adt_y(u_plr, rho_plr, lsp_y, flux_N, flux_S);

        // now we are ready to calculate the CU Numerical fluxes
        vvvd F = get_cu_x(u_plr, rho_plr, lsp_x, flux_E, flux_W, adt_x);
        vvvd G = get_cu_y(u_plr, rho_plr, lsp_y, flux_N, flux_S, adt_y);

        return {F,G};

    }

    // Helper functions to accomplish parts of the process are defined below:-

    // function to extend u and rho matrices on all boundaries
    void ExtendVarMatrices(){
        // extend cells
        EXC::extend_matrix(u, bc, "VELOCITY");
        EXC::extend_matrix(rho, bc, "DENSITY");
    }

    // Function to calculate new time step using the CFL conditions
    double update_dt(vvvd& lsp_x, vvvd& lsp_y){
        double amax=0.0; 
        double bmax=0.0;

        // both lsp_x and lsp_y have the same dimensions but different relevant values
        int row = lsp_x[0].size();
        int col = lsp_x[0][0].size();

        for(int i=1 ; i<row-1 ; i++){
            for(int j=0 ; j<col-1 ; j++){
                amax = max(
                    amax,
                    max(
                        lsp_x[0][i][j], -1*lsp_x[1][j][j]
                        // lsp_x[0][i][j], lsp_x[1][j][j]
                    )
                );
            }
        }


        for(int i=0 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                bmax = max(
                    bmax,
                    max(
                        lsp_y[0][i][j] , -1*lsp_y[1][i][j]
                        // lsp_y[0][i][j] , lsp_y[1][i][j]
                    )
                );
            }
        }

        // log
        cout << "amax (t=" << t << ") = " << amax << endl;
        cout << "bmax (t=" << t << ") = " << bmax << endl;

        double dt = CTS::CFL*min( dx/amax , dy/bmax );

        if( t+dt > time.second ){
            dt = time.second - t;
        }

        return dt;
    }

    // function to obtain the PLR of a given conserved variable
    // used as a subroutine in CU flux calculation
    // [N, S, E, W, NE, NW, SE, SW] -> respectively
    vvvd get_plr(vvd& var, string var_name){
        vvd slx = get_slx(var);        
        vvd sly = get_sly(var);        

        int row = var.size();
        int col = var[0].size();
        
        vvvd plr; // [N, S, E, W, NE, NW, SE, SW] -> respectively

        vvd n(row-2, vd(col-2));
        vvd s(row-2, vd(col-2));
        vvd e(row-2, vd(col-2));
        vvd w(row-2, vd(col-2));
        vvd ne(row-2, vd(col-2));
        vvd nw(row-2, vd(col-2));
        vvd se(row-2, vd(col-2));
        vvd sw(row-2, vd(col-2));

        for( int i=1 ; i<row-1 ; i++ ){
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

        // extend cells for the piecewise linear reconstructions
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

    // functions to get the slx and sly matrices
    // used as a subroutine in PLR calculation
    vvd get_slx(vvd& var){
        int row = var.size();
        int col = var[0].size();

        vvd slx(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                double v1 = CTS::THETA * (var[i][j] - var[i][j-1]);
                double v2 = 0.5*(var[i][j+1] - var[i][j-1]);
                double v3 = CTS::THETA * (var[i][j+1] - var[i][j]);

                slx[i-1][j-1] = UTL::minmod( v1 , v2 , v3 );
            }
        }

        EXC::extend_matrix(slx, bc, "SLOPES");

        return slx;
    }

    vvd get_sly(vvd& var){
        int row = var.size();
        int col = var[0].size();

        vvd sly(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                double v1 = CTS::THETA * (var[i][j] - var[i-1][j]);
                double v2 = 0.5*(var[i+1][j] - var[i-1][j]);
                double v3 = CTS::THETA * (var[i+1][j] - var[i][j]);

                sly[i-1][j-1] = UTL::minmod( v1 , v2 , v3 );
            }
        }

        EXC::extend_matrix(sly, bc, "SLOPES");

        return sly;
    }

    // function to calculate the One sided local speeds of propagations
    // a+, a-
    vvvd get_lsp_a(vvvd& u_plr, vvvd& rho_plr){
        int row = u_plr[0].size();
        int col = u_plr[0][0].size();

        vvvd lsp;

        vvd ap(row, vd(col));
        vvd am(row, vd(col));

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){

                ap[i][j] = max(
                    0.0 ,
                    max(
                        // abs( 2*u_plr[3][i][j+1] ),
                        // abs( 2*u_plr[2][i][j] )
                        2*u_plr[3][i][j+1],
                        2*u_plr[2][i][j]
                    )
                );

                am[i][j] = min(
                    0.0,
                    min(
                        // abs( u_plr[3][i][j+1] ),
                        // abs( u_plr[2][i][j] )
                        u_plr[3][i][j+1],
                        u_plr[2][i][j]
                    )
                );

            }
        }

        lsp.push_back(ap);
        lsp.push_back(am);

        // log
        // cout << "a+ (t=" << t << ")" << endl;
        // for(auto i : ap){
        //     for(auto j : i){
        //         cout << j << " ";
        //     }
        //     cout << endl;
        // }

        // cout << "a- (t=" << t << ")" << endl;
        // for(auto i : am){
        //     for(auto j : i){
        //         cout << j << " ";
        //     }
        //     cout << endl;
        // }

        return lsp;
    }   

    // function to calculate the One sided local speeds of propagations
    // b+, b-
    vvvd get_lsp_b(vvvd& u_plr, vvvd& rho_plr){
        int row = u_plr[0].size();
        int col = u_plr[0][0].size();

        vvvd lsp;

        vvd bp(row, vd(col));
        vvd bm(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){

                bp[i][j] = max(
                    0.0 ,
                    max(
                        // abs( 2*u_plr[1][i+1][j] ),
                        // abs( 2*u_plr[0][i][j] )
                        2*u_plr[1][i+1][j],
                        2*u_plr[0][i][j]
                    )
                );

                bm[i][j] = min(
                    0.0,
                    min(
                        // abs( u_plr[1][i+1][j] ),
                        // abs( u_plr[0][i][j] )
                        u_plr[1][i+1][j],
                        u_plr[0][i][j]
                    )
                );

            }
        }

        lsp.push_back(bp);
        lsp.push_back(bm);

        return lsp;
    }   

    // Functions to calculate the Anti-diffusion terms
    // u_adt, rho_adt
    vvvd get_adt_x(vvvd& u_plr, vvvd& rho_plr, vvvd& lsp, vvvd& flux_E, vvvd& flux_W){

        vvd u_star = get_star_x(u_plr, lsp, flux_E[0], flux_W[0]);
        vvd rho_star = get_star_x(rho_plr, lsp, flux_E[1], flux_W[1]);

        // calculating anti diffusion terms
        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();
        
        vvvd adt;

        vvd u_adt(row, vd(col));
        vvd rho_adt(row, vd(col));

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){
                u_adt[i][j] = UTL::minmod(
                    u_plr[7][i][j+1]-u_star[i][j],
                    u_star[i][j]-u_plr[6][i][j],
                    u_plr[5][i][j+1]-u_star[i][j],
                    u_star[i][j] - u_plr[4][i][j]
                );

                rho_adt[i][j] = UTL::minmod(
                    rho_plr[7][i][j+1]-rho_star[i][j],
                    rho_star[i][j]-rho_plr[6][i][j],
                    rho_plr[5][i][j+1]-rho_star[i][j],
                    rho_star[i][j] - rho_plr[4][i][j]
                );
            }
        }

        adt.push_back(u_adt);
        adt.push_back(rho_adt);

        return adt;
    }

    vvvd get_adt_y(vvvd& u_plr, vvvd& rho_plr, vvvd& lsp, vvvd& flux_N, vvvd& flux_S){

        vvd u_star = get_star_y(u_plr, lsp, flux_N[0], flux_S[0]);
        vvd rho_star = get_star_y(rho_plr, lsp, flux_N[1], flux_S[1]);

        // calculating anti diffusion terms
        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();
        
        vvvd adt;

        vvd rho_adt(row, vd(col));
        vvd u_adt(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){
                u_adt[i][j] = UTL::minmod(
                    u_plr[7][i+1][j]-u_star[i][j],
                    u_star[i][j]-u_plr[5][i][j],
                    u_plr[6][i+1][j]-u_star[i][j],
                    u_star[i][j] - u_plr[4][i][j]
                );

                rho_adt[i][j] = UTL::minmod(
                    rho_plr[7][i+1][j]-rho_star[i][j],
                    rho_star[i][j]-rho_plr[5][i][j],
                    rho_plr[6][i+1][j]-rho_star[i][j],
                    rho_star[i][j] - rho_plr[4][i][j]
                );

            }
        }

        adt.push_back(u_adt);
        adt.push_back(rho_adt);

        return adt;
    }

    // Functions to calculate the Intermediate star terms used in calculation for the Anti-diffusion terms
    vvd get_star_x(vvvd& plr, vvvd& lsp, vvd& flux_E, vvd& flux_W){
        int row = plr[0].size();
        int col = plr[0][0].size();

        vvd star(row, vd(col));

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){
                star[i][j] = (lsp[0][i][j]*plr[3][i][j+1] - lsp[1][i][j]*plr[2][i][j] - flux_W[i][j+1] + flux_E[i][j]) / (lsp[0][i][j]-lsp[1][i][j]);
            }
        }

        return star;
    }

    vvd get_star_y(vvvd& plr, vvvd& lsp, vvd& flux_N, vvd& flux_S){
        int row = plr[0].size();
        int col = plr[0][0].size();

        vvd star(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){
                star[i][j] = (lsp[0][i][j]*plr[1][i+1][j] - lsp[1][i][j]*plr[0][i][j] - flux_S[i+1][j] + flux_N[i][j]) / (lsp[0][i][j]-lsp[1][i][j]);
            }
        }

        return star;
    }

    // Functions to calculate the CU Numerical flux given all the intermediate values
    vvvd get_cu_x(vvvd& u_plr, vvvd& rho_plr, vvvd& lsp, vvvd& flux_E, vvvd& flux_W, vvvd& adt){
        // Calculate the CU Numerical flux using all the intermediate terms
        vvvd F;
        
        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();

        vvd F1(row, vd(col));
        vvd F2(row, vd(col));

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){

                double diff = lsp[0][i][j] - lsp[1][i][j];
                double prod = lsp[0][i][j]*lsp[1][i][j];
                
                if( diff > CTS::EPSILON ){
                    // use ADT
                    F1[i][j] = ( lsp[0][i][j]*flux_E[0][i][j] - lsp[1][i][j]*flux_W[0][i][j+1] + prod*( u_plr[3][i][j+1] - u_plr[2][i][j] - adt[0][i][j] ) ) / (diff);
                    F2[i][j] = ( lsp[0][i][j]*flux_E[1][i][j] - lsp[1][i][j]*flux_W[1][i][j+1] + prod*( rho_plr[3][i][j+1] - rho_plr[2][i][j] - adt[1][i][j] ) ) / (diff);
                }
                else{
                    // use normal averaging
                    F1[i][j] = 0.5*(flux_E[0][i][j]+flux_W[0][i][j]);
                    F2[i][j] = 0.5*(flux_E[1][i][j]+flux_W[1][i][j]);
                }

            }
        }

        F.push_back(F1);
        F.push_back(F2);

        return F;
    }

    vvvd get_cu_y(vvvd& u_plr, vvvd& rho_plr, vvvd& lsp, vvvd& flux_N, vvvd& flux_S, vvvd& adt){
        // Calculate the CU Numerical flux using all the intermediate terms
        vvvd G;
        
        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();

        vvd G1(row, vd(col));
        vvd G2(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){

                double diff = lsp[0][i][j] - lsp[1][i][j];
                double prod = lsp[0][i][j]*lsp[1][i][j];
                
                if( diff > CTS::EPSILON ){
                    // use ADT
                    G1[i][j] = ( lsp[0][i][j]*flux_N[0][i][j] - lsp[1][i][j]*flux_S[0][i+1][j] + prod*( u_plr[1][i+1][j] - u_plr[0][i][j] - adt[0][i][j] ) ) / (diff);
                    G2[i][j] = ( lsp[0][i][j]*flux_N[1][i][j] - lsp[1][i][j]*flux_S[1][i+1][j] + prod*( rho_plr[1][i+1][j] - rho_plr[0][i][j] - adt[1][i][j] ) ) / (diff);
                }
                else{
                    // use normal averaging
                    G1[i][j] = 0.5*(flux_N[0][i][j]+flux_S[0][i][j]);
                    G2[i][j] = 0.5*(flux_N[1][i][j]+flux_S[1][i][j]);

                }

            }
        }

        G.push_back(G1);
        G.push_back(G2);

        return G;
    }

};



