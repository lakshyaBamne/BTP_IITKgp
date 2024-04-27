#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<cmath>
#include<fstream>
#include<sstream>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define pss pair<string,string>
#define pdd pair<double,double>

namespace CTS{
    //! Global constants
    double GAMMA = 1.4;
    double THETA = 1.3;
    double EPS = 1.0E-12;

    //! CFL number
    // double CFL = 0.9; // 1st order schemes
    double CFL = 0.45; // 2nd order schemes

    //! domain endpoints
    pdd domainX = {0, 1};
    // pdd domainX = {-5, 5}; // LAX
    // pdd domainX = {-5, 15}; // SDW
    // pdd domainX = {-5,5}; // SEW

    //! time endpoints
    // pdd time = {0, 2}; // MCW
    pdd time = {0, 0.012}; // SCW
    // pdd time = {0, 0.038}; // BLW
    // pdd time = {0, 1.3}; // LAX
    // pdd time = {0, 0.2}; // TORO-1
    // pdd time = {0, 0.15}; // TORO-2
    // pdd time = {0, 0.012}; // TORO-3
    // pdd time = {0, 0.035}; // TORO-4
    // pdd time = {0, 5}; // SDW
    // pdd time = {0, 5}; // SEW

    //! problem 
    // string PROBLEM = "MCW";
    string PROBLEM = "SCW";
    // string PROBLEM = "BLW";
    // string PROBLEM = "LAX";
    // string PROBLEM = "TORO-1";
    // string PROBLEM = "TORO-2";
    // string PROBLEM = "TORO-3";
    // string PROBLEM = "TORO-4";
    // string PROBLEM = "SDW";
    // string PROBLEM = "SEW";

    string BC = "FREE";
    // string BC = "REFLECTIVE";

    //! number of finite volume cells
    int N = 1000;
    // int N = 10;
    int Nref = 4000;

    // length of an individual finite volume cell
    double dx = (domainX.second-domainX.first)/N;
}

//! Function definitions
vd make_grid(pdd domainX, int N);
vvd initialize_conserved_variables(string PROBLEM, string BC, vd& grid, vd& rho, vd& u, vd& p, vd& a);

// Function to update the conserved variables
void update_conserved_variables(vvd& U, vd& rho, vd& u, vd& p, vd& a, double& dt);

//! Utility functions used in the scheme
void extend_fvcells(string BC, vvd& U, vd& rho, vd& u, vd& p, vd& a);
void update_time_step(double& dt, vd& u, vd& a);
vvd get_slopes(vvd& U);
vvvd get_reconstruction(vvd& U, vvd& sigma);
void extend_reconstructions(string BC, vvd& Up, vvd& Um);

// Functions to calculate the reconstructed variables
vd get_density(vvd& U);
vd get_velocity(vvd& U);
vd get_pressure(vvd& U, vd& rho, vd& u);
vd get_acoustic_speed(vd& p, vd& rho);

// function to calculate the alpha term for LLF Fluxes
vd get_alpha(vd& up, vd& um, vd& ap, vd& am);
vvd get_llf_flux(vvd& flxp, vvd& flxm, vd& alpha, vvd& Up, vvd& Um);

vvd get_flux(vd& rho, vd& u, vd& p);

double minmod(double a, double b, double c);

// Functions to output the result
void write_grid(vd& grid, string file);
void write_density(vd& density, string file);

// Function to run the scheme
void run_llf_scheme(string MODE);

int main(){
    run_llf_scheme("NORMAL");
    // run_llf_scheme("REF");

    return 0;
}

//! Function bodies

// Function to initialize a computational grid
vd make_grid(pdd domainX, int N){
    double dx = (domainX.second-domainX.first)/N;

    vd grid(N+2);

    for(int i=1 ; i<=N ; i++){
        grid[i] = domainX.first + (i - 0.5)*dx;
    }

    return grid;
}

// Function to initialize the conserved variables
vvd initialize_conserved_variables(string PROBLEM, string BC, vd& grid, vd& rho, vd& u, vd& p, vd& a){
    
    int N = rho.size()-2;

    if(PROBLEM == "MCW"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.3){
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
    else if(PROBLEM == "SCW"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.8){
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
    else if(PROBLEM == "BLW"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.1){
                rho[i] = 1;
                u[i] = 0;
                p[i] = 1000;
            }
            else if(grid[i]>=0.1 && grid[i]<=0.9){
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
    else if(PROBLEM == "LAX"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0){
                rho[i] = 0.445;
                u[i] = 0.698;
                p[i] = 3.528;
            }
            else{
                rho[i] = 0.500;
                u[i] = 0;
                p[i] = 0.571;
            }
        }
    }
    else if(PROBLEM == "TORO-1"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.3){
                rho[i] = 1.0;
                u[i] = 0.75;
                p[i] = 1.0;
            }
            else{
                rho[i] = 0.125;
                u[i] = 0;
                p[i] = 0.1;
            }
        }
    }
    else if(PROBLEM == "TORO-2"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.5){
                rho[i] = 1.0;
                u[i] = -2.0;
                p[i] = 0.4;
            }
            else{
                rho[i] = 1.0;
                u[i] = 2.0;
                p[i] = 0.4;
            }
        }
    }
    else if(PROBLEM == "TORO-3"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.5){
                rho[i] = 1.0;
                u[i] = 0;
                p[i] = 1000;
            }
            else{
                rho[i] = 1.0;
                u[i] = 0;
                p[i] = 0.01;
            }
        }
    }
    else if(PROBLEM == "TORO-4"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.4){
                rho[i] = 5.99924;
                u[i] = 19.5975;
                p[i] = 460.894;
            }
            else{
                rho[i] = 5.99242;
                u[i] = -6.19633;
                p[i] = 46.0950;
            }
        }
    }
    else if(PROBLEM == "SDW"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < -4){
                rho[i] = 27/4;
                u[i] = (4/9)*sqrt(35);
                p[i] = 31/3;
            }
            else{
                rho[i] = 1+0.2*sin(5*grid[i]);
                u[i] = 0;
                p[i] = 1;
            }
        }
    }
    else if(PROBLEM == "SEW"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < -4.5){
                rho[i] = 1.51695;
                u[i] = 0.523346;
                p[i] = 1.805;
            }
            else{
                rho[i] = 1+0.1*sin(20*grid[i]);
                u[i] = 0;
                p[i] = 1;
            }
        }
    }
    else{
        cout << "---LOG--- Please enter correct problem---" << endl;
    }

    // now initialize the conserved variables
    vvd U(3, vd(N+2, 0));

    for(int i=1 ; i<=N ; i++){
        U[0][i] = rho[i];
        U[1][i] = rho[i]*u[i];
        U[2][i] = (p[i]/(CTS::GAMMA-1)) + (0.5*rho[i]*pow(u[i],2));

        a[i] = sqrt(CTS::GAMMA*p[i]/rho[i]);
    }

    return U;
}

// Function to write the computational grid to external file
void write_grid(vd& grid, string file){
    stringstream ss;
    for(int i=1 ; i<grid.size()-1 ; i++){
        ss << grid[i] << " ";
    }
    
    ofstream fout;
    fout.open(file, ios::app);

    fout << ss.str() << endl;

    fout.close();
}

// Function to write the density to an external file
void write_density(vd& density, string file){
    stringstream ss;
    for(int i=1 ; i<density.size()-1 ; i++){
        ss << density[i] << " ";
    }

    ofstream fout;
    fout.open(file, ios::app);

    fout << ss.str() << endl;

    fout.close();
}

// Function to update conserved variables
void update_conserved_variables(vvd& U, vd& rho, vd& u, vd& p, vd& a, double& dt){

    // extend the cells for later stages
    extend_fvcells(CTS::BC, U, rho, u, p, a);

    //! update the time step
    update_time_step(dt, u, a);

    // get the reconstruction slopes
    vvd sigma = get_slopes(U);

    // get the 2nd order reconstructions for the conserved variables
    vvvd recons = get_reconstruction(U, sigma);
    vvd Up = recons[0];
    vvd Um = recons[1];

    // get the required fluxes using the reconstructions
    vd rho_p = get_density(Up);
    vd rho_m = get_density(Um);

    vd u_p = get_velocity(Up);
    vd u_m = get_velocity(Um);

    vd p_p = get_pressure(Up, rho_p, u_p);
    vd p_m = get_pressure(Um, rho_m, u_m);

    vd a_p = get_acoustic_speed(p_p, rho_p);
    vd a_m = get_acoustic_speed(p_m, rho_m);

    vvd flux_p = get_flux(rho_p, u_p, p_p);
    vvd flux_m = get_flux(rho_m, u_m, p_m);

    //! Now we are ready to calculate the alpha term and finally the LLF Fluxes
    vd alpha = get_alpha(u_p, u_m, a_p, a_m);
    vvd llf_flux = get_llf_flux(flux_p, flux_m, alpha, Up, Um);

    // update the conserved variables
    for(int i=0 ; i<3 ; i++){
        for(int j=1 ; j<=CTS::N ; j++){
            U[i][j] = U[i][j] - (dt/CTS::dx)*(llf_flux[i][j] - llf_flux[i][j-1]);
        }
    }

}

// Function to extend the ghost cells
void extend_fvcells(string BC, vvd& U, vd& rho, vd& u, vd& p, vd& a){
    if(BC=="FREE"){
        // extend conserved variables
        for(int i=0 ; i<3 ; i++){
            U[i][0] = U[i][1];
            U[i][CTS::N+1] = U[i][CTS::N];
        }

        // extend other variables
        rho[0] = rho[1];
        rho[CTS::N+1] = rho[CTS::N];

        u[0] = u[1];
        u[CTS::N+1] = u[CTS::N];

        p[0] = p[1];
        p[CTS::N+1] = p[CTS::N];

        a[0] = a[1];
        a[CTS::N+1] = a[CTS::N];

    }
    else{
        cout << "---ERROR--- Please select correct boundary conditions---" << endl;
    }
}

void extend_reconstructions(string BC, vvd& Up, vvd& Um){
    if(BC=="FREE"){
        for(int i=0 ; i<3 ; i++){
            Up[i][0] = Um[i][1];
            Um[i][CTS::N] = Up[i][CTS::N-1];
        }
    }
    else{
        cout << "---ERROR--- Please select correct boundary conditions---" << endl;
    }
}

// Function to update the time step using the CFL conditions and local speeds of propagation
void update_time_step(double& dt, vd& u, vd& a){
    double amax;

    for(int i=1 ; i<=CTS::N ; i++){
        amax = max(
            amax,
            abs(u[i]) + a[i]
        );
    }

    dt = CTS::CFL*CTS::dx/amax;
}

double minmod(double a, double b, double c){
    double max_val = max(a, max(b,c));
    double min_val = min(a, min(b,c));

    if(max_val < 0){
        return max_val;
    }
    else if(min_val > 0){
        return min_val;
    }
    else{
        return 0;
    }
}

// Get the slopes for the internal grid points
vvd get_slopes(vvd& U){
    vvd sigma(3, vd(CTS::N+2));

    for(int i=0 ; i<3 ; i++){
        for(int j=1 ; j<=CTS::N ; j++){
            sigma[i][j] = minmod(
                CTS::THETA*(U[i][j] - U[i][j-1]),
                0.5*(U[i][j+1] - U[i][j-1]),
                CTS::THETA*(U[i][j+1] - U[i][j])
            );
        }
    }

    return sigma;
}

// Get the polynomial reconstruction
vvvd get_reconstruction(vvd& U, vvd& sigma){
    vvd Up(3, vd(CTS::N+1));
    vvd Um(3, vd(CTS::N+1));

    for(int i=0 ; i<3 ; i++){
        // U+
        for(int j=1 ; j<=CTS::N ; j++){
            Up[i][j] = U[i][j] + 0.5*sigma[i][j];
        }

        // U-
        for(int j=0 ; j<=CTS::N-1 ; j++){
            Um[i][j] = U[i][j+1] - 0.5*sigma[i][j+1];
        }
    }

    // extend the reconstructions
    extend_reconstructions(CTS::BC, Up, Um);

    return {Up, Um};
}

vvd get_flux(vd& rho, vd& u, vd& p){
    vvd F(3, vd(CTS::N+1));

    for(int i=0 ; i<=CTS::N ; i++){
        F[0][i] = rho[i]*u[i];
        F[1][i] = rho[i]*pow(u[i],2) + p[i];
        F[2][i] = u[i]*(p[i] + ((p[i]/(CTS::GAMMA-1)) + (0.5*rho[i]*pow(u[i],2))));
    }

    return F;
}

vd get_velocity(vvd& U){
    vd u(CTS::N+1);

    for(int i=0 ; i<=CTS::N ; i++){
        u[i] = U[1][i]/U[0][i];
    }

    return u;
}

vd get_pressure(vvd& U, vd& rho, vd& u){
    vd p(CTS::N+1);

    for(int i=0 ; i<=CTS::N ; i++){
        p[i] = (CTS::GAMMA-1)*(U[2][i] - 0.5*(rho[i]*pow(u[i],2)));
    }

    return p;
}

vd get_acoustic_speed(vd& p, vd& rho){
    vd a(CTS::N+1);

    for(int i=0 ; i<=CTS::N ; i++){
        a[i] = sqrt(CTS::GAMMA*p[i]/rho[i]);
    }

    return a;
}

vd get_density(vvd& U){
    vd rho(CTS::N+1);

    for(int i=0 ; i<=CTS::N ; i++){
        rho[i] = U[0][i];
    }

    return rho;
}

vd get_alpha(vd& up, vd& um, vd& ap, vd& am){
    vd alpha(CTS::N+1);

    for(int i=0 ; i<=CTS::N ; i++){
        alpha[i] = max(
            abs(um[i])+am[i],
            abs(up[i])+ap[i]
        );
    }

    return alpha;
}

vvd get_llf_flux(vvd& flxp, vvd& flxm, vd& alpha, vvd& Up, vvd& Um){
    vvd llf(3, vd(CTS::N+1));

    for(int i=0 ; i<3 ; i++){
        for(int j=0 ; j<=CTS::N ; j++){
            llf[i][j] = 0.5*(flxp[i][j] + flxm[i][j]) - 0.5*alpha[j]*(Up[i][j] - Um[i][j]);
        }
    }

    return llf;
}

//! Function to run the scheme
void run_llf_scheme(string MODE){
    if(MODE=="NORMAL"){
        //! Initialize computational grid
        vd grid = make_grid(CTS::domainX, CTS::N);

        //! write the computational grid to an output text file
        write_grid(grid, MODE+"grid.txt");

        //! initialize the conserved variables vector based on the given problem
        vd rho(CTS::N+2), u(CTS::N+2), p(CTS::N+2), a(CTS::N+2);
        vvd U = initialize_conserved_variables(CTS::PROBLEM, CTS::BC, grid, rho, u, p, a);

        double dt;
        double t=CTS::time.first;

        while(t < CTS::time.second){
            // log message
            cout << "t=" << t << "|dt=" << dt << endl;

            //! write to output file
            write_density(U[0], MODE+"density.txt");

            update_conserved_variables(U, rho, u, p, a, dt);
            t += dt;
        }

    }
    else{
        cout << "---ERROR--- Please select a mode for running the scheme---" << endl;
    }
}

