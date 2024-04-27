#include<iostream>

#include "Constants.h"
#include "OutputResult.h"

using namespace std;

namespace cts = Constants;
namespace opr = OutputResult;

double minmod(double a, double b);
void initialize_conserved_variables(string problem, vd& grid, vd& U1, vd& U2, vd& U3);
void extend_boundaries(string bc, vd& U1, vd& U2, vd& U3);
void update_conserved_variables(vd& U1, vd& U2, vd& U3, double& dt);

int main(){
    int n = cts::n;

    //! other variables and arrays used in the simulations
    double dt, t;

    vd x(n+2); // computational grid

    vd U1(n+2); // density
    vd U2(n+2); // momentum
    vd U3(n+2); // energy

    //! Initialize the computational grid
    for(int i=1 ; i<=n ; i++){
        x[i] = cts::domx.first + cts::dx*(i-0.5);
    }

    // output computational grid
    opr::write_vector(x, "grid.txt");

    //! Initialie the Riemann Problem (primitve variables)
    initialize_conserved_variables(cts::problem, x, U1, U2, U3);

    //! extend boundaries for the primitive variables
    extend_boundaries(cts::bc, U1, U2, U3);

    // output initial state
    opr::write_vector(U1, "density.txt");

    t = cts::time.first;
    while( t < cts::time.second ){
        cout << "t=" << t << " | dt=" << dt << endl;
        update_conserved_variables(U1, U2, U3, dt);

        t += dt;
    }

    // output final state
    opr::write_vector(U1, "density.txt");

    cout << "------------------------------" << endl;
    cout << "Final Time = " << t << endl;
    cout << "------------------------------" << endl;

    return 0;
}

void update_conserved_variables(vd& U1, vd& U2, vd& U3, double& dt){

    int n = U1.size() - 2;

    //! Calculate the primitive variables using the Conserved variables for use later
    vd density(n+2);
    vd velocity(n+2);
    vd pressure(n+2);
    vd acoustic_speed(n+2); // speed of sound

    for(int i=0 ; i<=n+1 ; i++){
        density[i] = U1[i];
        velocity[i] = U2[i]/U1[i];
        pressure[i] = (cts::GAMMA-1)*( U3[i] - pow(U2[i],2)/(2*U1[i]) );
        acoustic_speed[i] = sqrtf( cts::GAMMA*pressure[i]/density[i] );
    }

    // arguments for the minmod limiter
    vd sigma1(n+2);
    vd sigma2(n+2);
    vd sigma3(n+2);

    //! Left of the interface
    // slopes to the left of the (j+1/2)th interface
    vd slx_left1(n+2);
    vd slx_left2(n+2);
    vd slx_left3(n+2);

    // second order reconstructions left of (j+1/2)th interface
    vd Uleft1(n+2);
    vd Uleft2(n+2);
    vd Uleft3(n+2);

    //! Right of the interface
    // slopes to the right of the (j+1/2)th interface
    vd slx_right1(n+2);
    vd slx_right2(n+2);
    vd slx_right3(n+2);

    // second order reconstructions right of (j+1/2)th interface
    vd Uright1(n+2);
    vd Uright2(n+2);
    vd Uright3(n+2);

    // flux vector at left of interface
    vd fleft1(n+2);
    vd fleft2(n+2);
    vd fleft3(n+2);

    // flux vector at right of interface
    vd fright1(n+2);
    vd fright2(n+2);
    vd fright3(n+2);

    // intermediate terms used in the LLF Flux
    vd alpha(n+2);

    // Vectors to store the LLF Fluxes at the (j+1/2)th interface and (j-1/2)th interface
    vd llf_ph1(n+2);
    vd llf_ph2(n+2);
    vd llf_ph3(n+2);

    vd llf_mh1(n+2);
    vd llf_mh2(n+2);
    vd llf_mh3(n+2);

    /*
        ! LLF Flux for the (j+1/2)th interface
    */

    // arguments for the minmod limiter calculated for each grid point (using Forward Differences here)
    for(int i=1 ; i<=n ; i++){
        sigma1[i] = U1[i+1] - U1[i];
        sigma2[i] = U2[i+1] - U2[i];
        sigma3[i] = U3[i+1] - U3[i];
    }

    //! Calculate slopes left of the j+1/2 interface present at the jth grid point (3, n-2)
    for(int i=3 ; i<=n-2 ; i++){
        slx_left1[i] = minmod(sigma1[i-1], sigma1[i]);
        slx_left2[i] = minmod(sigma2[i-1], sigma2[i]);
        slx_left3[i] = minmod(sigma3[i-1], sigma3[i]);
    }

    for(int i=3 ; i<=n-2 ; i++){
        Uleft1[i] = U1[i] + 0.5*slx_left1[i];
        Uleft2[i] = U2[i] + 0.5*slx_left2[i];
        Uleft3[i] = U3[i] + 0.5*slx_left3[i];
    }

    //! Calculate slopes right of the j+1/2 interface present at the jth grid point (3, n-2)
    for(int i=3 ; i<=n-2 ; i++){
        slx_right1[i] = minmod(sigma1[i+1], sigma1[i]);
        slx_right2[i] = minmod(sigma2[i+1], sigma2[i]);
        slx_right3[i] = minmod(sigma3[i+1], sigma3[i]);
    }

    for(int i=3 ; i<=n-2 ; i++){
        Uright1[i] = U1[i+1] - 0.5*slx_right1[i];
        Uright2[i] = U2[i+1] - 0.5*slx_right2[i];
        Uright3[i] = U3[i+1] - 0.5*slx_right3[i];
    }

    // fluxes at left of (j+1/2)th interface
    for(int i=3 ; i<=n-2 ; i++){
        fleft1[i] = Uleft2[i];
        fleft2[i] = ( pow(Uleft2[i],2)/Uleft1[i] ) + (cts::GAMMA-1)*( Uleft3[i] - ( pow(Uleft2[i],2)/(2*Uleft1[i]) ) );
        fleft3[i] = ( Uleft2[i]/Uleft1[i] )*( (cts::GAMMA*Uleft3[i]) - (cts::GAMMA-1)*( pow(Uleft2[i],2)/(2*Uleft1[i]) ) );
    }

    // fluxes at right of (j+1/2)th interface
    for(int i=3 ; i<=n-2 ; i++){
        fright1[i] = Uright2[i];
        fright2[i] = ( pow(Uright2[i],2)/Uright1[i] ) + (cts::GAMMA-1)*( Uright3[i] - ( pow(Uright2[i],2)/(2*Uright1[i]) ) );
        fright3[i] = ( Uright2[i]/Uright1[i] )*( (cts::GAMMA*Uright3[i]) - (cts::GAMMA-1)*( pow(Uright2[i],2)/(2*Uright1[i]) ) );
    }

    // intermediate alpha term calculation for left and right of the (j+1/2)th interface
    for(int i=3 ; i<=n ; i++){
        alpha[i] = max(
            abs(velocity[i]) + acoustic_speed[i],
            abs(velocity[i+1]) + acoustic_speed[i+1]
        );
    }

    // LLF FLux for the (j+1/2)th interface
    for(int i=3 ; i<=n-2 ; i++){
        llf_ph1[i] = 0.5*( fright1[i] + fleft1[i] ) - 0.5*alpha[i]*( Uright1[i] - Uleft1[i] );
        llf_ph2[i] = 0.5*( fright2[i] + fleft2[i] ) - 0.5*alpha[i]*( Uright2[i] - Uleft2[i] );
        llf_ph3[i] = 0.5*( fright3[i] + fleft3[i] ) - 0.5*alpha[i]*( Uright3[i] - Uleft3[i] );
    }

    /*
        ! LLF Flux for the (j-1/2)th interface
    */

    // arguments for the minmod limiter calculated for each grid point (using Forward Differences here)
    for(int i=1 ; i<=n ; i++){
        sigma1[i] = U1[i] - U1[i-1];
        sigma2[i] = U2[i] - U2[i-1];
        sigma3[i] = U3[i] - U3[i-1];
    }

    //! Calculate slopes left of the j-1/2 interface present at the jth grid point (3, n-2)
    for(int i=3 ; i<=n-2 ; i++){
        slx_left1[i] = minmod(sigma1[i-1], sigma1[i]);
        slx_left2[i] = minmod(sigma2[i-1], sigma2[i]);
        slx_left3[i] = minmod(sigma3[i-1], sigma3[i]);
    }

    for(int i=3 ; i<=n-2 ; i++){
        Uleft1[i] = U1[i-1] + 0.5*slx_left1[i];
        Uleft2[i] = U2[i-1] + 0.5*slx_left2[i];
        Uleft3[i] = U3[i-1] + 0.5*slx_left3[i];
    }

    //! Calculate slopes right of the j-1/2 interface present at the jth grid point (3, n-2)
    for(int i=3 ; i<=n-2 ; i++){
        slx_right1[i] = minmod(sigma1[i+1], sigma1[i]);
        slx_right2[i] = minmod(sigma2[i+1], sigma2[i]);
        slx_right3[i] = minmod(sigma3[i+1], sigma3[i]);
    }

    for(int i=3 ; i<=n-2 ; i++){
        Uright1[i] = U1[i] - 0.5*slx_right1[i];
        Uright2[i] = U2[i] - 0.5*slx_right2[i];
        Uright3[i] = U3[i] - 0.5*slx_right3[i];
    }

    // fluxes at left of (j-1/2)th interface
    for(int i=3 ; i<=n-2 ; i++){
        fleft1[i] = Uleft2[i];
        fleft2[i] = ( pow(Uleft2[i],2)/Uleft1[i] ) + (cts::GAMMA-1)*( Uleft3[i] - ( pow(Uleft2[i],2)/(2*Uleft1[i]) ) );
        fleft3[i] = ( Uleft2[i]/Uleft1[i] )*( (cts::GAMMA*Uleft3[i]) - (cts::GAMMA-1)*( pow(Uleft2[i],2)/(2*Uleft1[i]) ) );
    }

    // fluxes at right of (j-1/2)th interface
    for(int i=3 ; i<=n-2 ; i++){
        fright1[i] = Uright2[i];
        fright2[i] = ( pow(Uright2[i],2)/Uright1[i] ) + (cts::GAMMA-1)*( Uright3[i] - ( pow(Uright2[i],2)/(2*Uright1[i]) ) );
        fright3[i] = ( Uright2[i]/Uright1[i] )*( (cts::GAMMA*Uright3[i]) - (cts::GAMMA-1)*( pow(Uright2[i],2)/(2*Uright1[i]) ) );
    }

    // intermediate alpha term calculation for left and right of the (j-1/2)th interface
    for(int i=3 ; i<=n ; i++){
        alpha[i] = max(
            abs(velocity[i-1]) + acoustic_speed[i-1],
            abs(velocity[i]) + acoustic_speed[i]
        );
    }

    // LLF FLux for the (j-1/2)th interface
    for(int i=3 ; i<=n-2 ; i++){
        llf_mh1[i] = 0.5*( fright1[i] + fleft1[i] ) - 0.5*alpha[i]*( Uright1[i] - Uleft1[i] );
        llf_mh2[i] = 0.5*( fright2[i] + fleft2[i] ) - 0.5*alpha[i]*( Uright2[i] - Uleft2[i] );
        llf_mh3[i] = 0.5*( fright3[i] + fleft3[i] ) - 0.5*alpha[i]*( Uright3[i] - Uleft3[i] );
    }


    /*
        ! Update the time step
    */
    double amax = 0;
    for(int i=1 ; i<=n ; i++){
        amax = max(amax, abs(velocity[i])+acoustic_speed[i]);
    }
    dt = cts::CFL*cts::dx/amax;

    /*
        ! Update the conserved variables using the LLF Fluxes for next time step using Euler Forward Differences
    */
    vd U1new(n+2);
    vd U2new(n+2);
    vd U3new(n+2);

    for(int i=3 ; i<=n-2 ; i++){
        U1new[i] = U1[i] - (dt/cts::dx)*( llf_ph1[i] - llf_mh1[i] );
        U2new[i] = U2[i] - (dt/cts::dx)*( llf_ph2[i] - llf_mh2[i] );
        U3new[i] = U3[i] - (dt/cts::dx)*( llf_ph3[i] - llf_mh3[i] );
    }

    // update conserved variables
    for(int i=3 ; i<=n-2 ; i++){
        U1[i] = U1new[i];
        U2[i] = U2new[i];
        U3[i] = U3new[i];
    }


}

double minmod(double a, double b){
    if( min(a,b) > 0 ){
        return min(a,b);
    }        
    else if( max(a,b) < 0 ){
        return max(a,b);
    }
    else{
        return 0;
    }
}

void initialize_conserved_variables(string problem, vd& grid, vd& U1, vd& U2, vd& U3){
    
    int n = grid.size() - 2;
    double u, p;
    
    if( problem == "SCW" ){

        for(int i=1 ; i<=n ; i++){
            
            U1[i] = 1.0;
            u = -19.59745;

            if( grid[i] <= 0.8 ){
                p = 1000.0;
            }
            else{
                p = 0.01;
            }

            U2[i] = U1[i]*u;
            U3[i] = (p/(cts::GAMMA-1)) + 0.5*U1[i]*pow(u,2);
        }

    }
    else{
        cout << "---ERROR--- Please enter correct problem---" << endl;
    }

}

void extend_boundaries(string bc, vd& U1, vd& U2, vd& U3){
    int n = U1.size()-2;
    
    if( bc == "FREE" ){
        // density
        U1[0] = U1[1];
        U1[n+1] = U1[n];
        // momentum
        U2[0] = U2[1];
        U2[n+1] = U2[n];
        // energy
        U3[0] = U3[1];
        U3[n+1] = U3[n];

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions---" << endl;
    }
}