#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>

#include "Constants.h"
#include "ExtendCells.h"
#include "Utility.h"
#include "OutputResult.h"
#include "InitializeConservedVariables.h"

using namespace std;

namespace cts = Constants;
namespace exc = ExtendCells;
namespace utl = Utility;
namespace opr = OutputResult;
namespace icv = InitializeConservedVariables;

//! Function definitions
void export_env();
vvvd get_plr(vvd U, string var);
vvvvd get_cu_flux(vvvd U, bool update_time, double& dt, double T);
void update_conserved_variables(vvvd& U, double& dt, double T);

int main(){
    /*
        ! Export the run environment for later use
    */
    export_env();

    /*
        ! Definitions for the required variables
    */
    double T = cts::time.first;
    double dt;

    //! Define computaitonal grids and the conserved variables vector
    // computational grid vectors
    vd gridx(cts::nx+2,0);
    vd gridy(cts::ny+2,0);

    // conserved variable vectors
    vvvd U(4, vvd(cts::nx+2, vd(cts::ny+2, 0))); // [rho, mx, my, e]

    /*
        ! Simulation steps start from now
    */

    //! Initialize the computational grid lines
    for(int i=1 ; i<=cts::nx ; i++){
        gridx[i] = cts::domx.first + cts::dx*( i - 0.5 );
    }

    for(int j=1 ; j<=cts::ny ; j++){
        gridy[j] = cts::domy.first + cts::dy*( j - 0.5 );
    }

    // output the computational grids
    opr::write_grids(cts::MODE, gridx, gridy);

    /*
        ! Initialize the conserved variables and then extend the cells for ghost values
    */
    icv::initialize_conserved_variables(U, gridx, gridy);
    exc::extend_conserved_variables(U);

    // output initial conditions
    opr::write_matrix(cts::MODE, "Density", U[0]); // write density
    
    /*
        ! Now we can start the iterations
    */
    while( T < cts::time.second ){
        cout << "T=" << T << " | dt=" << dt << endl;
        update_conserved_variables(U, dt, T);
        T += dt;
    }

    // output final conditions
    opr::write_matrix(cts::MODE, "Density", U[0]); // write density

    return 0;

}

//! Function implementations

/*
    Export the run environment for later use in plotting
*/
void export_env(){
    string file_name;

    if( cts::MODE == 1 ){
        file_name = "env1/Environment.txt";
    }
    else{
        file_name = "env2/Environment.txt";
    }
    
    stringstream ss;
    
    ofstream fout;
    fout.open(file_name , ios::app);

    ss << cts::MODE << endl;
    ss << cts::PROBLEM << endl;
    ss << cts::domx.first << " " << cts::domx.second << endl;
    ss << cts::domy.first << " " << cts::domy.second << endl;
    ss << cts::time.first << " " << cts::time.second << endl;
    ss << cts::nx << endl;
    ss << cts::ny << endl;

    fout << ss.str() << endl;

    fout.close();
}

/*
    Update conserved variables
*/
void update_conserved_variables(vvvd& U, double& dt, double T){
    
    //! Stage-1
    vvvvd flux = get_cu_flux(U, true, dt, T);
    vvvd F = flux[0];
    vvvd G = flux[1];

    // get the parameters used in updation
    double LAMBDA = dt / cts::dx;
    double MU = dt / cts::dy;

    vvvd U1(4, vvd(cts::nx+2, vd(cts::ny+2, 0))); // [rho, mx, my, e]
    
    for(int i=1 ; i<=cts::nx ; i++){
        for(int j=1 ; j<=cts::ny ; j++){

            for(int u=0 ; u<4 ; u++){
                U1[u][i][j] = U[u][i][j] - ( LAMBDA*( F[u][i][j] - F[u][i-1][j] ) + MU*( G[u][i][j] - G[u][i][j-1] ) );
            }

        }
    }

    exc::extend_conserved_variables(U1);

    //! Stage-2
    vvvvd flux2 = get_cu_flux(U1, false, dt, T);
    vvvd F1 = flux2[0];
    vvvd G1 = flux2[1];
    
    vvvd U2(4, vvd(cts::nx+2, vd(cts::ny+2, 0))); // [rho, mx, my, e]
    
    for(int i=1 ; i<=cts::nx ; i++){
        for(int j=1 ; j<=cts::ny ; j++){

            for(int u=0 ; u<4 ; u++){
                U2[u][i][j] = ( (3.0*U[u][i][j] + U1[u][i][j]) - ( LAMBDA*( F1[u][i][j] - F1[u][i-1][j] ) + MU*( G1[u][i][j] - G1[u][i][j-1] ) ) ) / 4.0;
            }

        }
    }

    exc::extend_conserved_variables(U2);

    //! Stage-3
    vvvvd flux3 = get_cu_flux(U2, false, dt, T);
    vvvd F2 = flux3[0];
    vvvd G2 = flux3[1];
    
    vvvd UN(4, vvd(cts::nx+2, vd(cts::ny+2, 0))); // [rho, mx, my, e]
    
    for(int i=1 ; i<=cts::nx ; i++){
        for(int j=1 ; j<=cts::ny ; j++){

            for(int u=0 ; u<4 ; u++){
                UN[u][i][j] = ( ( U[u][i][j] + 2.0*U2[u][i][j] ) - 2.0*( LAMBDA*( F2[u][i][j] - F2[u][i-1][j] ) + MU*( G2[u][i][j] - G2[u][i][j-1] ) ) ) / 3.0;
            }

        }
    }

    exc::extend_conserved_variables(UN);

    // update the original vector containing the conserved variables
    // U = UN;
    for(int u=0 ; u<4 ; u++){
        for(int i=0 ; i<cts::nx+2 ; i++){
            for(int j=0 ; j<cts::ny+2 ; j++){
                U[u][i][j] = U1[u][i][j];
            }
        }
    }

}

/*
    Function to calculate the CU Numerical Flux given a set of conserved variables
*/
vvvvd get_cu_flux(vvvd U, bool update_time, double& dt, double T){
    /*
        ! First we need to get the Piecewise Linear Reconstructions
    */
    vvvvd PLR;

    PLR.push_back( get_plr(U[0], "Density") );
    PLR.push_back( get_plr(U[1], "MomentumX") );
    PLR.push_back( get_plr(U[2], "MomentumY") );
    PLR.push_back( get_plr(U[3], "Energy") );

    /*
        ! Now we can calculate the CU Numerical Flux directly with 
        ! local speeds of propagation as an intermediate
    */

    double uE, uW, vE, vW, pE, pW;
    double f1E, f2E, f3E, f4E;
    double f1W, f2W, f3W, f4W;

    // horizontal direction
    vvd F1(cts::nx+2, vd(cts::ny+2, 0));
    vvd F2(cts::nx+2, vd(cts::ny+2, 0));
    vvd F3(cts::nx+2, vd(cts::ny+2, 0));
    vvd F4(cts::nx+2, vd(cts::ny+2, 0));

    vvd ap(cts::nx+2, vd(cts::ny+2, 0));
    vvd am(cts::nx+2, vd(cts::ny+2, 0));
    
    for(int i=0 ; i<=cts::nx ; i++){
        for(int j=1 ; j<=cts::ny ; j++){
            uE = PLR[1][2][i][j] / PLR[0][2][i][j];
            uW = PLR[1][3][i+1][j] / PLR[0][3][i+1][j];

            vE = PLR[2][2][i][j] / PLR[0][2][i][j];
            vW = PLR[2][3][i+1][j] / PLR[0][3][i+1][j];
            
            pE = ( cts::GAMMA - 1 )*( PLR[3][2][i][j] - 0.5*PLR[0][2][i][j]*( pow(uE,2) + pow(vE,2) ) );
            pW = ( cts::GAMMA - 1 )*( PLR[3][3][i+1][j] - 0.5*PLR[0][3][i+1][j]*( pow(uW,2) + pow(vW,2) ) );

            ap[i][j] = max(
                0.0,
                max(
                    uW + sqrt( cts::GAMMA*pW / PLR[0][3][i+1][j] ),
                    uE + sqrt( cts::GAMMA*pE / PLR[0][2][i][j] )
                )
            );

            am[i][j] = min(
                0.0,
                min(
                    uW - sqrt( cts::GAMMA*pW / PLR[0][3][i+1][j] ),
                    uE - sqrt( cts::GAMMA*pE / PLR[0][2][i][j] )
                )
            );

            f1E = PLR[1][2][i][j];
            f2E = PLR[1][2][i][j]*uE + pE;
            f3E = PLR[1][2][i][j]*vE;
            f4E = uE*( PLR[3][2][i][j] + pE );

            f1W = PLR[1][3][i+1][j];
            f2W = PLR[1][3][i+1][j]*uW + pW;
            f3W = PLR[1][3][i+1][j]*vW;
            f4W = uW*( PLR[3][3][i+1][j] + pW );

            double dist = ap[i][j] - am[i][j];
            double prod = ap[i][j]*am[i][j];

            if( dist > cts::EPSILON ){
                F1[i][j] = ( ( ap[i][j]*f1E - am[i][j]*f1W ) + ( prod * ( PLR[0][3][i+1][j] - PLR[0][2][i][j] ) ) ) / dist;
                F2[i][j] = ( ( ap[i][j]*f2E - am[i][j]*f2W ) + ( prod * ( PLR[1][3][i+1][j] - PLR[1][2][i][j] ) ) ) / dist;
                F3[i][j] = ( ( ap[i][j]*f3E - am[i][j]*f3W ) + ( prod * ( PLR[2][3][i+1][j] - PLR[2][2][i][j] ) ) ) / dist;
                F4[i][j] = ( ( ap[i][j]*f4E - am[i][j]*f4W ) + ( prod * ( PLR[3][3][i+1][j] - PLR[3][2][i][j] ) ) ) / dist;
            }
            else{
                F1[i][j] = 0.5 * ( f1E + f1W );
                F2[i][j] = 0.5 * ( f2E + f2W );
                F3[i][j] = 0.5 * ( f3E + f3W );
                F4[i][j] = 0.5 * ( f4E + f4W );
            }

        }
    }

    // vertical direction
    double uN, uS, vN, vS, pN, pS;
    double g1N, g2N, g3N, g4N;
    double g1S, g2S, g3S, g4S;

    vvd G1(cts::nx+2, vd(cts::ny+2, 0));
    vvd G2(cts::nx+2, vd(cts::ny+2, 0));
    vvd G3(cts::nx+2, vd(cts::ny+2, 0));
    vvd G4(cts::nx+2, vd(cts::ny+2, 0));

    vvd bp(cts::nx+2, vd(cts::ny+2, 0));
    vvd bm(cts::nx+2, vd(cts::ny+2, 0));

    for(int i=1 ; i<=cts::nx ; i++){
        for(int j=0 ; j<=cts::ny ; j++){

            uN = PLR[1][0][i][j] / PLR[0][0][i][j];
            uS = PLR[1][1][i][j+1] / PLR[0][1][i][j+1];
            vN = PLR[2][0][i][j] / PLR[0][0][i][j];
            vS = PLR[2][1][i][j+1] / PLR[0][1][i][j+1];
            pN = ( cts::GAMMA - 1 )*( PLR[3][0][i][j] - 0.5*PLR[0][0][i][j]*( pow(uN,2) + pow(vN,2) ) );
            pS = ( cts::GAMMA - 1 )*( PLR[3][1][i][j+1] - 0.5*PLR[0][1][i][j+1]*( pow(uS,2) + pow(vS,2) ) );

            bp[i][j] = max(
                0.0,
                max(
                    vS + sqrt( cts::GAMMA*pS / PLR[0][1][i][j+1] ),
                    vN + sqrt( cts::GAMMA*pN / PLR[0][0][i][j] )
                )
            );

            bm[i][j] = min(
                0.0,
                min(
                    vS - sqrt( cts::GAMMA*pS / PLR[0][1][i][j+1] ),
                    vN - sqrt( cts::GAMMA*pN / PLR[0][0][i][j] )
                )
            );

            g1N = PLR[2][0][i][j];
            g2N = PLR[2][0][i][j]*uN;
            g3N = PLR[2][0][i][j]*vN + pN;
            g4N = vN*( PLR[3][0][i][j] + pN );

            g1S = PLR[2][1][i][j+1];
            g2S = PLR[2][1][i][j+1]*uS;
            g3S = PLR[2][1][i][j+1]*vS + pS;
            g4S = vS*( PLR[3][1][i][j+1] + pS );

            double dist = bp[i][j] - bm[i][j];
            double prod = bp[i][j] * bm[i][j];

            if( dist > cts::EPSILON ){
                G1[i][j] = ( ( bp[i][j]*g1N - bm[i][j]*g1S ) + ( prod * ( PLR[0][1][i][j+1] - PLR[0][1][i][j] ) ) ) / dist;
                G2[i][j] = ( ( bp[i][j]*g2N - bm[i][j]*g2S ) + ( prod * ( PLR[1][1][i][j+1] - PLR[1][1][i][j] ) ) ) / dist;
                G3[i][j] = ( ( bp[i][j]*g3N - bm[i][j]*g3S ) + ( prod * ( PLR[2][1][i][j+1] - PLR[2][1][i][j] ) ) ) / dist;
                G4[i][j] = ( ( bp[i][j]*g4N - bm[i][j]*g4S ) + ( prod * ( PLR[3][1][i][j+1] - PLR[3][1][i][j] ) ) ) / dist;
            }
            else{
                G1[i][j] = 0.5 * ( g1N + g1S );
                G2[i][j] = 0.5 * ( g2N + g2S );
                G3[i][j] = 0.5 * ( g3N + g3S );
                G4[i][j] = 0.5 * ( g4N + g4S );
            }

        }
    }

    /*
        ! Now we can update the time step using the CFL conditions and the local speeds of propagation
    */
    if( update_time ){
        
        double amax = 0.0;
        double bmax = 0.0;

        for(int i=0 ; i<=cts::nx ; i++){
            for(int j=1 ; j<=cts::ny ; j++){
                amax = max(
                    amax,
                    max(ap[i][j], -am[i][j])
                );
            }
        }

        for(int i=1 ; i<=cts::nx ; i++){
            for(int j=0 ; j<=cts::ny ; j++){
                bmax = max(
                    bmax,
                    max(bp[i][j], -bm[i][j])
                );
            }
        }

        dt = cts::CFL * min( cts::dx/amax , cts::dy/bmax );

        if( T+dt > cts::time.second ){
            dt = cts::time.second - T;
        }

    }

    return {{F1, F2, F3, F4}, {G1, G2, G3, G4}};

}

/*
    Function to calculate the Piecewise Linear Reconstructions
*/
vvvd get_plr(vvd U, string var){

    vvd n(cts::nx+2, vd(cts::ny+2, 0)); // 0
    vvd s(cts::nx+2, vd(cts::ny+2, 0)); // 1
    vvd e(cts::nx+2, vd(cts::ny+2, 0)); // 2
    vvd w(cts::nx+2, vd(cts::ny+2, 0)); // 3

    double slx, sly;

    for(int i=1 ; i<=cts::nx ; i++){
        for(int j=1 ; j<=cts::ny ; j++){
            slx = utl::minmod( cts::THETA*( U[i][j] - U[i-1][j] ) , 0.5*(U[i+1][j] - U[i-1][j]) , cts::THETA*( U[i+1][j] - U[i][j] ) );
            sly = utl::minmod( cts::THETA*( U[i][j] - U[i][j-1] ) , 0.5*(U[i][j+1] - U[i][j-1]) , cts::THETA*( U[i][j+1] - U[i][j] ) );

            n[i][j] = U[i][j] + 0.5*sly;
            s[i][j] = U[i][j] - 0.5*sly;
            e[i][j] = U[i][j] + 0.5*slx;
            w[i][j] = U[i][j] - 0.5*slx;
        }
    }

    // extend the plr
    exc::extend_plr(n, s, e, w, var);

    return {n, s, e, w};

}

