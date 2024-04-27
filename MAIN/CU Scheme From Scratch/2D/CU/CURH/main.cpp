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
vvvd get_plr(vvd U);
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
        update_conserved_variables(U, dt, T);
        T += dt;
        cout << "T=" << T << " | dt=" << dt << endl;
    }

    // output final conditions
    opr::write_matrix(cts::MODE, "Density", U[0]); // write density

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
    
    for(int j=1 ; j<=cts::ny ; j++){
        for(int i=1 ; i<=cts::nx ; i++){

            for(int u=0 ; u<4 ; u++){
                U1[u][i][j] = U[u][i][j] - ( LAMBDA*( F[u][i][j] - F[u][i-1][j] ) + MU*( G[u][i][j] - G[u][i][j-1] ) );
            }

        }
    }

    exc::extend_conserved_variables(U1);

    //! Stage-2
    flux = get_cu_flux(U1, false, dt, T);
    F = flux[0];
    G = flux[1];
    
    vvvd U2(4, vvd(cts::nx+2, vd(cts::ny+2, 0))); // [rho, mx, my, e]
    
    for(int j=1 ; j<=cts::ny ; j++){
        for(int i=1 ; i<=cts::nx ; i++){

            for(int u=0 ; u<4 ; u++){
                U2[u][i][j] = ( (3.0*U[u][i][j] + U1[u][i][j]) - ( LAMBDA*( F[u][i][j] - F[u][i-1][j] ) + MU*( G[u][i][j] - G[u][i][j-1] ) ) ) / 4.0;
            }

        }
    }

    exc::extend_conserved_variables(U2);

    //! Stage-3
    flux = get_cu_flux(U2, false, dt, T);
    F = flux[0];
    G = flux[1];
    
    vvvd UN(4, vvd(cts::nx+2, vd(cts::ny+2, 0))); // [rho, mx, my, e]
    
    for(int j=1 ; j<=cts::ny ; j++){
        for(int i=1 ; i<=cts::nx ; i++){

            for(int u=0 ; u<4 ; u++){
                UN[u][i][j] = ( ( U[u][i][j] + 2.0*U2[u][i][j] ) - 2.0*( LAMBDA*( F[u][i][j] - F[u][i-1][j] ) + MU*( G[u][i][j] - G[u][i][j-1] ) ) ) / 3.0;
            }

        }
    }

    exc::extend_conserved_variables(UN);

    // update the original vector containing the conserved variables
    U = UN;

}

/*
    Function to calculate the CU Numerical Flux given a set of conserved variables
*/
vvvvd get_cu_flux(vvvd U, bool update_time, double& dt, double T){
    /*
        ! First we need to get the Piecewise Linear Reconstructions
    */
    vvvvd PLR;
    
    for(int u=0 ; u<4 ; u++){
        PLR.push_back( get_plr(U[u]) );
    }

    /*
        ! Now we can calculate the CU Numerical Flux directly with 
        ! local speeds of propagation as an intermediate
    */

    // horizontal direction
    vvd F1(cts::nx+2, vd(cts::ny+2, 0));
    vvd F2(cts::nx+2, vd(cts::ny+2, 0));
    vvd F3(cts::nx+2, vd(cts::ny+2, 0));
    vvd F4(cts::nx+2, vd(cts::ny+2, 0));

    vvd ap(cts::nx+2, vd(cts::ny+2, 0));
    vvd am(cts::nx+2, vd(cts::ny+2, 0));
    
    for(int j=1 ; j<=cts::ny ; j++){
        for(int i=0 ; i<=cts::nx ; i++){
            double uE = PLR[1][2][i][j] / PLR[0][2][i][j];
            double uW = PLR[1][3][i+1][j] / PLR[0][3][i+1][j];
            double vE = PLR[2][2][i][j] / PLR[0][2][i][j];
            double vW = PLR[2][3][i+1][j] / PLR[0][3][i+1][j];
            double pE = ( cts::GAMMA - 1 )*( PLR[3][2][i][j] - 0.5*PLR[0][2][i][j]*( pow(uE,2) + pow(vE,2) ) );
            double pW = ( cts::GAMMA - 1 )*( PLR[3][3][i+1][j] - 0.5*PLR[0][3][i+1][j]*( pow(uW,2) + pow(vW,2) ) );

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

            /*
                ! CURH Optimization block-1 added by Prof. Naveen Garg
            */
            if( cts::ScaleSw ){
                double differ1 = abs( 0.5*(PLR[0][3][i+1][j]*pow(uW,2) - PLR[0][2][i][j]*pow(uE,2)) + ((pW-pE)/(cts::GAMMA-1)) );
                double differ2 = 0.5*abs(PLR[0][3][i+1][j]*pow(vW,2) - PLR[0][2][i][j]*pow(vE,2));
                double differ = sqrt( pow(differ1,2) + pow(differ2,2) );

                double alpha;

                if( differ > cts::EPSILON ){
                    alpha = differ1/differ;
                }
                else{
                    alpha = 0;
                }

                ap[i][j] = max(
                    0.0,
                    max(
                        uW + alpha*sqrt( cts::GAMMA*pW/PLR[0][3][i+1][j] ),
                        uE + alpha*sqrt( cts::GAMMA*pE/PLR[0][2][i][j] )
                    )
                );

                ap[i][j] = min(
                    0.0,
                    min(
                        uW - alpha*sqrt( cts::GAMMA*pW/PLR[0][3][i+1][j] ),
                        uE - alpha*sqrt( cts::GAMMA*pE/PLR[0][2][i][j] )
                    )
                );
            }

            double f1E = PLR[1][2][i][j];
            double f2E = PLR[1][2][i][j]*uE + pE;
            double f3E = PLR[1][2][i][j]*vE;
            double f4E = uE*( PLR[3][2][i][j] + pE );

            double f1W = PLR[1][3][i+1][j];
            double f2W = PLR[1][3][i+1][j]*uW + pW;
            double f3W = PLR[1][3][i+1][j]*vW;
            double f4W = uW*( PLR[3][3][i+1][j] + pW );

            /*
                ! CURH Optimization block-2 added by Prof. Naveen Garg
            */
            if( cts::ShrinkSw ){
                double df1 = f1W - f1E;
                double df2 = f2W - f2E;
                double df3 = f3W - f3E;
                double df4 = f4W - f4E;
            
                double du1 = PLR[0][3][i+1][j] - PLR[0][2][i][j];
                double du2 = PLR[1][3][i+1][j] - PLR[1][2][i][j];
                double du3 = PLR[2][3][i+1][j] - PLR[2][2][i][j];
                double du4 = PLR[3][3][i+1][j] - PLR[3][2][i][j];

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

                double ratio1 = ( 2.0*df1 ) / ( du1 + du1eps );
                double ratio2 = ( 2.0*df2 ) / ( du2 + du2eps );
                double ratio3 = ( 2.0*df3 ) / ( du3 + du3eps );
                double ratio4 = ( 2.0*df4 ) / ( du4 + du4eps );

                double ratiomin = min( ratio1, min(ratio2, min(ratio3, ratio4)) );
                double ratiomax = max( ratio1, max(ratio2, max(ratio3, ratio4)) );

                if( ratiomax > 0 ){
                    ap[i][j] = min( ap[i][j] , ratiomax );
                    am[i][j] = max( am[i][j] , -ratiomax );
                }

                if( ratiomin < 0 ){
                    ap[i][j] = min( ap[i][j] , -ratiomin );
                    am[i][j] = max( am[i][j] , ratiomin );
                }

            }

            double dist = ap[i][j] - am[i][j];
            double prod = ap[i][j]*am[i][j];

            if( dist > cts::EPSILON ){
                /*
                    ! Calculating Anti-diffusion term for better solution evolution
                    ! -> this code block was included to form the new version of the CU Scheme
                    * Requirement is to calculate the Piecewise Linear Approxmiations for NE, NW, SE, SW also
                */
                double rho_star = ( ( ap[i][j]*PLR[0][3][i+1][j] - am[i][j]*PLR[0][1][i][j] ) - ( f1W - f1E ) ) / dist;
                double mx_star = ( ( ap[i][j]*PLR[1][3][i+1][j] - am[i][j]*PLR[1][1][i][j] ) - ( f2W - f2E ) ) / dist;
                double my_star = ( ( ap[i][j]*PLR[2][3][i+1][j] - am[i][j]*PLR[2][1][i][j] ) - ( f3W - f3E ) ) / dist;
                double E_star = ( ( ap[i][j]*PLR[3][3][i+1][j] - am[i][j]*PLR[3][1][i][j] ) - ( f4W - f4E ) ) / dist;

                double drho = utl::minmod(
                    PLR[0][7][i+1][j] - rho_star,
                    rho_star - PLR[0][6][i][j],
                    PLR[0][5][i+1][j] - rho_star,
                    rho_star - PLR[0][4][i][j]
                );

                double dmx = utl::minmod(
                    PLR[1][7][i+1][j] - mx_star,
                    mx_star - PLR[1][6][i][j],
                    PLR[1][5][i+1][j] - mx_star,
                    mx_star - PLR[1][4][i][j]
                );

                double dmy = utl::minmod(
                    PLR[2][7][i+1][j] - my_star,
                    my_star - PLR[2][6][i][j],
                    PLR[2][5][i+1][j] - my_star,
                    my_star - PLR[2][4][i][j]
                );

                double dE = utl::minmod(
                    PLR[3][7][i+1][j] - E_star,
                    E_star - PLR[3][6][i][j],
                    PLR[3][5][i+1][j] - E_star,
                    E_star - PLR[3][4][i][j]
                );

                F1[i][j] = ( ( ap[i][j]*f1E - am[i][j]*f1W ) + ( prod * ( PLR[0][3][i+1][j] - PLR[0][2][i][j] - drho ) ) ) / dist;
                F2[i][j] = ( ( ap[i][j]*f2E - am[i][j]*f2W ) + ( prod * ( PLR[1][3][i+1][j] - PLR[1][2][i][j] - dmx ) ) ) / dist;
                F3[i][j] = ( ( ap[i][j]*f3E - am[i][j]*f3W ) + ( prod * ( PLR[2][3][i+1][j] - PLR[2][2][i][j] - dmy ) ) ) / dist;
                F4[i][j] = ( ( ap[i][j]*f4E - am[i][j]*f4W ) + ( prod * ( PLR[3][3][i+1][j] - PLR[3][2][i][j] - dE ) ) ) / dist;
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
    vvd G1(cts::nx+2, vd(cts::ny+2, 0));
    vvd G2(cts::nx+2, vd(cts::ny+2, 0));
    vvd G3(cts::nx+2, vd(cts::ny+2, 0));
    vvd G4(cts::nx+2, vd(cts::ny+2, 0));

    vvd bp(cts::nx+2, vd(cts::ny+2, 0));
    vvd bm(cts::nx+2, vd(cts::ny+2, 0));

    for(int j=0 ; j<=cts::ny ; j++){
        for(int i=1 ; i<=cts::nx ; i++){

            double uN = PLR[1][0][i][j] / PLR[0][0][i][j];
            double uS = PLR[1][1][i][j+1] / PLR[0][1][i][j+1];
            double vN = PLR[2][0][i][j] / PLR[0][0][i][j];
            double vS = PLR[2][1][i][j+1] / PLR[0][1][i][j+1];
            double pN = ( cts::GAMMA - 1 )*( PLR[3][0][i][j] - 0.5*PLR[0][0][i][j]*( pow(uN,2) + pow(vN,2) ) );
            double pS = ( cts::GAMMA - 1 )*( PLR[3][1][i][j+1] - 0.5*PLR[0][1][i][j+1]*( pow(uS,2) + pow(vS,2) ) );

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

            /*
                ! CURH Optimization block-1 added by Prof. Naveen Garg
            */
            if( cts::ScaleSw ){
                double differ1 = abs( 0.5*(PLR[0][1][i][j+1]*pow(vS,2) - PLR[0][0][i][j]*pow(vN,2)) + ((pS-pN)/(cts::GAMMA-1)) );
                double differ2 = 0.5*abs(PLR[0][1][i][j+1]*pow(uS,2) - PLR[0][0][i][j]*pow(uN,2));
                double differ = sqrt( pow(differ1,2) + pow(differ2,2) );

                double alpha;

                if( differ > cts::EPSILON ){
                    alpha = differ1/differ;
                }
                else{
                    alpha = 0;
                }

                bp[i][j] = max(
                    0.0,
                    max(
                        vS + alpha*sqrt( cts::GAMMA*pS/PLR[0][1][i][j+1] ),
                        vN + alpha*sqrt( cts::GAMMA*pN/PLR[0][0][i][j] )
                    )
                );

                bm[i][j] = min(
                    0.0,
                    min(
                        vS - alpha*sqrt( cts::GAMMA*pS/PLR[0][1][i][j+1] ),
                        vN - alpha*sqrt( cts::GAMMA*pN/PLR[0][0][i][j] )
                    )
                );
            }

            double g1N = PLR[2][0][i][j];
            double g2N = PLR[2][0][i][j]*uN;
            double g3N = PLR[2][0][i][j]*vN + pN;
            double g4N = vN*( PLR[3][0][i][j] + pN );

            double g1S = PLR[2][1][i][j+1];
            double g2S = PLR[2][1][i][j+1]*uS;
            double g3S = PLR[2][1][i][j+1]*vS + pS;
            double g4S = vS*( PLR[3][1][i][j+1] + pS );

            /*
                ! CURH Optimization block-2 added by Prof. Naveen Garg
            */
            if( cts::ShrinkSw ){
                double dg1 = g1S - g1N;
                double dg2 = g2S - g2N;
                double dg3 = g3S - g3N;
                double dg4 = g4S - g4N;
            
                double du1 = PLR[0][1][i][j+1] - PLR[0][0][i][j];
                double du2 = PLR[1][1][i][j+1] - PLR[1][0][i][j];
                double du3 = PLR[2][1][i][j+1] - PLR[2][0][i][j];
                double du4 = PLR[3][1][i][j+1] - PLR[3][0][i][j];

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

                double ratio1 = ( 2.0*dg1 ) / ( du1 + du1eps );
                double ratio2 = ( 2.0*dg2 ) / ( du2 + du2eps );
                double ratio3 = ( 2.0*dg3 ) / ( du3 + du3eps );
                double ratio4 = ( 2.0*dg4 ) / ( du4 + du4eps );

                double ratiomin = min( ratio1, min(ratio2, min(ratio3, ratio4)) );
                double ratiomax = max( ratio1, max(ratio2, max(ratio3, ratio4)) );

                if( ratiomax > 0 ){
                    bp[i][j] = min( bp[i][j] , ratiomax );
                    bm[i][j] = max( bm[i][j] , -ratiomax );
                }

                if( ratiomin < 0 ){
                    bp[i][j] = min( bp[i][j] , -ratiomin );
                    bm[i][j] = max( bm[i][j] , ratiomin );
                }

            }

            double dist = bp[i][j] - bm[i][j];
            double prod = bp[i][j]*bm[i][j];

            if( dist > cts::EPSILON ){
                /*
                    ! Calculating Anti-diffusion term for better solution evolution
                    ! -> this code block was included to form the new version of the CU Scheme
                    * Requirement is to calculate the Piecewise Linear Approxmiations for NE, NW, SE, SW also
                */
                double rho_star = ( ( bp[i][j]*PLR[0][1][i][j+1] - bm[i][j]*PLR[0][0][i][j] ) - ( g1S - g1N ) ) / dist;
                double mx_star = ( ( bp[i][j]*PLR[1][1][i][j+1] - bm[i][j]*PLR[1][0][i][j] ) - ( g2S - g2N ) ) / dist;
                double my_star = ( ( bp[i][j]*PLR[2][1][i][j+1] - bm[i][j]*PLR[2][0][i][j] ) - ( g3S - g3N ) ) / dist;
                double E_star = ( ( bp[i][j]*PLR[3][1][i][j+1] - bm[i][j]*PLR[3][0][i][j] ) - ( g4S - g4N ) ) / dist;

                double drho = utl::minmod(
                    PLR[0][7][i][j+1] - rho_star,
                    rho_star - PLR[0][5][i][j],
                    PLR[0][6][i][j+1] - rho_star,
                    rho_star - PLR[0][4][i][j]
                );

                double dmx = utl::minmod(
                    PLR[1][7][i][j+1] - mx_star,
                    mx_star - PLR[1][5][i][j],
                    PLR[1][6][i][j+1] - mx_star,
                    mx_star - PLR[1][4][i][j]
                );

                double dmy = utl::minmod(
                    PLR[2][7][i][j+1] - my_star,
                    my_star - PLR[2][5][i][j],
                    PLR[2][6][i][j+1] - my_star,
                    my_star - PLR[2][4][i][j]
                );

                double dE = utl::minmod(
                    PLR[3][7][i][j+1] - E_star,
                    E_star - PLR[3][5][i][j],
                    PLR[3][6][i][j+1] - E_star,
                    E_star - PLR[3][4][i][j]
                );

                G1[i][j] = ( ( bp[i][j]*g1N - bm[i][j]*g1S ) + ( prod * ( PLR[0][1][i][j+1] - PLR[0][1][i][j] - drho ) ) ) / dist;
                G2[i][j] = ( ( bp[i][j]*g2N - bm[i][j]*g2S ) + ( prod * ( PLR[1][1][i][j+1] - PLR[1][1][i][j] - dmx ) ) ) / dist;
                G3[i][j] = ( ( bp[i][j]*g3N - bm[i][j]*g3S ) + ( prod * ( PLR[2][1][i][j+1] - PLR[2][1][i][j] - dmy ) ) ) / dist;
                G4[i][j] = ( ( bp[i][j]*g4N - bm[i][j]*g4S ) + ( prod * ( PLR[3][1][i][j+1] - PLR[3][1][i][j] - dE ) ) ) / dist;
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
        
        double amax = 0;
        double bmax = 0;

        for(int j=1 ; j<=cts::ny ; j++){
            for(int i=0 ; i<=cts::nx ; i++){
                amax = max(
                    amax,
                    max(ap[i][j], -am[i][j])
                );
            }
        }

        for(int j=0 ; j<=cts::ny ; j++){
            for(int i=1 ; i<=cts::nx ; i++){
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
vvvd get_plr(vvd U){

    vvvd plr;

    vvd n(cts::nx+2, vd(cts::ny+2, 0)); // 0
    vvd s(cts::nx+2, vd(cts::ny+2, 0)); // 1
    vvd e(cts::nx+2, vd(cts::ny+2, 0)); // 2
    vvd w(cts::nx+2, vd(cts::ny+2, 0)); // 3
    vvd ne(cts::nx+2, vd(cts::ny+2, 0)); // 4
    vvd nw(cts::nx+2, vd(cts::ny+2, 0)); // 5
    vvd se(cts::nx+2, vd(cts::ny+2, 0)); // 6
    vvd sw(cts::nx+2, vd(cts::ny+2, 0)); // 7

    double v1, v2, v3;
    double slx, sly;

    for(int i=1 ; i<=cts::nx ; i++){
        for(int j=1 ; j<=cts::ny ; j++){
            slx = utl::minmod( cts::THETA*( U[i][j] - U[i-1][j] ) , 0.5*(U[i+1][j] - U[i-1][j]) , cts::THETA*( U[i+1][j] - U[i][j] ) );
            sly = utl::minmod( cts::THETA*( U[i][j] - U[i][j-1] ) , 0.5*(U[i][j+1] - U[i][j-1]) , cts::THETA*( U[i][j+1] - U[i][j] ) );

            n[i][j] = U[i][j] + 0.5*sly;
            s[i][j] = U[i][j] - 0.5*sly;
            e[i][j] = U[i][j] + 0.5*slx;
            w[i][j] = U[i][j] - 0.5*slx;
            ne[i][j] = U[i][j] + 0.5*slx + 0.5*sly;
            nw[i][j] = U[i][j] - 0.5*slx + 0.5*sly;
            se[i][j] = U[i][j] + 0.5*slx - 0.5*sly;
            sw[i][j] = U[i][j] - 0.5*slx - 0.5*sly;
        }
    }

    // extend the plr
    exc::extend_plr(n, s, e, w, ne, nw, se, sw);

    return {n, s, e, w, ne, nw, se, sw};

}

