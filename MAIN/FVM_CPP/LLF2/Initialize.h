#pragma once

#include "Constants.h"

namespace CTS = Constants;

namespace Initialize{
    // Function to initialize the computational grid
    vd make_grid();   

    // Function to initialize the primitive variables
    vvd set_primitive_variables(const vd& x);

    // Function to extend the primitive variables
    void extend_primitive_variables(vd& rho, vd& u, vd& p);

    // Function to initialize the conserved variables 
    void initialize_conserved_variables(const vvd& P, vvd& U, vd& a);
}

vd Initialize::make_grid(){
    vd grid(CTS::N+2);

    for(int i=1 ; i<=CTS::N ; i++){
        grid[i] = CTS::domainX.first + (i - 0.5)*CTS::dx;
    }

    return grid;
}

vvd Initialize::set_primitive_variables(const vd& x){
    vd rho(CTS::N+2);
    vd u(CTS::N+2);
    vd p(CTS::N+2);

    if(CTS::PROBLEM == "MCW"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(x[i] < 0.3){
                rho[i] = 1.4;
                u[i] = 0.1;
                p[i] = 1;
            }
            else{
                rho[i] = 1.0;
                u[i] = 0.1;
                p[i] = 1;
            }
        }
    }
    else if(CTS::PROBLEM == "SCW"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(x[i] < 0.8){
                rho[i] = 1;
                u[i] = -19.59745;
                p[i] = 1000;
            }
            else{
                rho[i] = 1;
                u[i] = -19.59745;
                p[i] = 0.01;
            }
        }
    }
    else if(CTS::PROBLEM == "BLW"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(x[i] < 0.1){
                rho[i] = 1;
                u[i] = 0;
                p[i] = 1000;
            }
            else if(x[i]>=0.1 && x[i]<=0.9){
                rho[i] = 1;
                u[i] = 0;
                p[i] = 0.01;
            }
            else{
                rho[i] = 1;
                u[i] = 0;
                p[i] = 100;
            }
        }
    }
    else{
        cout << "---ERROR--- Please select correct problem---" << endl;
    }

    // extend the primitive variables
    extend_primitive_variables(rho, u, p);

    return {rho, u, p};
}

void Initialize::extend_primitive_variables(vd& rho, vd& u, vd& p){
    if(CTS::BC == "FREE"){
        rho[0] = rho[1];
        rho[CTS::N+1] = rho[CTS::N];

        u[0] = u[1];
        u[CTS::N+1] = u[CTS::N];

        p[0] = p[1];
        p[CTS::N+1] = p[CTS::N];
    }
    else if(CTS::BC == "REFLECTIVE"){
        rho[0] = rho[1];
        rho[CTS::N+1] = rho[CTS::N];

        u[0] = -u[1];
        u[CTS::N+1] = -u[CTS::N];

        p[0] = p[1];
        p[CTS::N+1] = p[CTS::N];
    }
    else{
        cout << "---ERROR--- Please select correct boundary conditions---" << endl;
    }
}

void Initialize::initialize_conserved_variables(const vvd& P, vvd& U, vd& a){
    for(int i=0 ; i<=CTS::N+1 ; i++){
        // conserved variables
        U[0][i] = P[0][i];
        U[1][i] = P[0][i]*P[1][i];
        U[2][i] = (P[2][i]/(CTS::GAMMA-1)) + 0.5*P[0][i]*pow(P[1][i],2);
    
        // acoustic speed
        a[i] = sqrt(CTS::GAMMA*U[2][i]/U[0][i]);
    }
}
