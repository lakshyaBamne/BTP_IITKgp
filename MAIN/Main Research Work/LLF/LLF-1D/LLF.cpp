/*
    ! Local Lax Friedrichs Scheme
*/

#include<iostream>
#include<cmath>

#include "Constants.h"

using namespace std;
namespace cts = Constants;

vd make_grid();
void initialize_variables(vd& x, vd& U1, vd& U2, vd& U3, vd& density, vd& velocity, vd& pressure, vd& acoustic_speed);
void update_time(vd& acoustic_speed, vd& velocity, double& dt);\
void extend_boundaries(vd& U1, vd& U2, vd& U3, vd& density, vd& velocity, vd& pressure, vd& acoustic_speed);

int main(){
    //! STEP-1 Initialize the computational grid
    vd x = make_grid();

    //! STEP-2 Initialize the conserved variables & extend boundaries
    vd density(cts::n+2, 0), velocity(cts::n+2, 0), pressure(cts::n+2, 0), acoustic_speed(cts::n+2, 0);
    vd U1(cts::n+2, 0), U2(cts::n+2, 0), U3(cts::n+2, 0);

    initialize_variables(x, U1, U2, U3, density, velocity, pressure, acoustic_speed);
    extend_boundaries(U1, U2, U3, density, velocity, pressure, acoustic_speed);

    double t, dt;
    update_time(acoustic_speed, velocity, dt);

    return 0;
}

// function to initialize the computational grid
vd make_grid(){
    vd grid(cts::n+2, 0);

    for(int i=1 ; i<=cts::n ; i++){
        grid[i] = cts::domx.first + cts::dx*( i-0.5 );
    }

    return grid;
}

// function to initialize the conserved variables
void initialize_variables(vd& x, vd& U1, vd& U2, vd& U3, vd& density, vd& velocity, vd& pressure, vd& acoustic_speed){

    if( cts::problem == "SCW" ){

        for(int i=1 ; i<=cts::n ; i++){
            
            // initialize the three primitive variables
            if( x[i] <= 0.8 ){
                density[i] = 1.0;
                velocity[i] = -19.59745;
                pressure[i] = 1000.0;
            }
            else{
                density[i] = 1.0;
                velocity[i] = -19.59745;
                pressure[i] = 0.01;
            }

            // initialize the conserved variables
            U1[i] = density[i];
            U2[i] = density[i]*velocity[i];
            U3[i] = (pressure[i]/(cts::GAMMA-1)) + (density[i]*pow(velocity[i],2)/2);

            // acoustic speed
            acoustic_speed[i] = sqrt(cts::GAMMA*pressure[i]/density[i]);
        }

    }
    else{
        cout << "---ERROR--- Please select correct problem ---" << endl;
    }

}

void update_time(vd& acoustic_speed, vd& velocity, double& dt){
    double amax = 0.0;
    for(int i=1 ; i<=cts::n ; i++){
        amax = max(amax, abs(velocity[i])+acoustic_speed[i]);
    }
    dt = cts::CFL*cts::dx/amax;
}

void extend_boundaries(vd& U1, vd& U2, vd& U3, vd& density, vd& velocity, vd& pressure, vd& acoustic_speed){
    if( cts::bc == "FREE" ){
        //! left boundary
        U1[0] = U1[1];
        U2[0] = U2[1];
        U3[0] = U3[1];

        density[0] = density[1];
        velocity[0] = velocity[1];
        pressure[0] = pressure[1];

        acoustic_speed[0] = acoustic_speed[1];

        //! right boundary
        U1[cts::n+1] = U1[cts::n];
        U2[cts::n+1] = U2[cts::n];
        U3[cts::n+1] = U3[cts::n];

        density[cts::n+1] = density[cts::n];
        velocity[cts::n+1] = velocity[cts::n];
        pressure[cts::n+1] = pressure[cts::n];

        acoustic_speed[cts::n+1] = acoustic_speed[cts::n];
    }
    else{
        cout << "---ERROR--- Please select correct boundary conditions ---" << endl;
    }
}
