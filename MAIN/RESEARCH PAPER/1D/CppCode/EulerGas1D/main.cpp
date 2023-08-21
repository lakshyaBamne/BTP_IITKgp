#include<iostream>
#include<vector>
#include<utility>
#include<algorithm>
#include<cmath>

#include "Constants.h" // Namespace containing all the constants being used in the program
#include "Utility.h" // Namespace containing useful Utility functions
#include "PiecewiseConstruction.h" // Namespace to calculate the Piecewise Linear Reconstruction
#include "SpeedsPropagation.h" // Namespace to calculate the One sided local speeds of propagation
#include "Flux.h" // Namespace to calculate the Flux vectors using the Conserved variables
#include "AntiDiffusionCalculator.h" // Namespace to calculate the Built in Anti-Diffusion term
#include "CuNumericalFLux.h" // Namespace to calculate the CU Numerical Flux
#include "CFL.h" // Namespace to calculate the new dt for the next iteration using the CFL conditions
#include "TimeDiscretization.h" // Namespace containing the 3-stage SSPRK scheme

#define ll long long int

using namespace std;

// namespace aliases
namespace CTS = Constants;
namespace PLR = PiecewiseLinearReconstruction;
namespace LSP = LocalSpeedsOfPropagation;
namespace FLX = Flux;
namespace ADC = AntiDuffusionCalculator;
namespace CUF = CUNumericalFlux;
namespace CFL = CFL;
namespace SRK = ThreeStageSSPRK;

int main(){
    /*
        Variable initializations required in the program

        -> variable definitions
        -> user input
    */
    
    double x_start; // start of computational domain
    double x_end; // end of computational domain
    
    ll Nx; // number of grid points in the computational domain
    double x_len; // length of the computational domain

    double dx; // individual cell length in the computational domain

    double T_final; // final time for iterations
    double T_initial; // initial time from where iterations start
    double dT; // time steps to go from initial to final time during iterations

    cout << "-----------------INPUT INITIAL CONDITIONS--------------" << endl;
    
    cout << "Domain start : ";
    cin >> x_start;

    cout << "Domain end : ";
    cin >> x_end;

    cout << "Grid points in domain : ";
    cin >> Nx;

    cout << "Initial time : ";
    cin >> T_initial;

    cout << "Final time : ";
    cin >> T_final;

    // computation of some variables based on the input given by the user
    
    // length of an individual cell in the Finite Volume grid
    x_len = x_end-x_start;
    dx = x_len / Nx;

    /*
        Data structures needed in the program are initialized here
    */

    vector<double> x(Nx); // computational grid

    vector<double> rho(Nx), u(Nx), p(Nx); // Primitive variables - density, velocity, proessure
    vector<double> m(Nx), E(Nx); // Conserved variables - momentum, energy

    //! Algorithm starts
    /*
        STEP-1 
        -> 1.1 Initialize the Finite volume grid for the computational domain
        -> 1.2 Initialize the Vectors containing values for the Conserved variables
    
    */
    cout << "-------------------------------------------------------" << endl;    
    cout << "-------------INITIALIZING FINITE VOLUME GRID-----------" << endl;
    cout << "-------------------------------------------------------" << endl;

    // 1.1 Initializing the Finite volume grid for the computational domain
    //! Information is stored in the geometric mid-point of a Finite Volume Grid Cell
    //! and flux refers to the quantities moving to and fro from the edges of these cells
    for(int i=1 ; i<=Nx ; i++ ){
        auto val = x_start + (i-0.5)*dx;
        x[i-1] = val;
    }

    // Log
    cout << "Finite Volume Grid (len=" << x.size() << ")" << endl;
    for(auto i : x){
        cout << i << " ";
    }
    cout << endl;

    cout << "-------------------------------------------------------" << endl;    
    cout << "-------------INITIALIZING CONSERVED VARIABLES----------" << endl;
    cout << "-------------------------------------------------------" << endl;

    // 1.2 Initializing the vectors containing values for the conserved variables
    
    // This can be changed based on the given Riemann problem
    //! Slowly moving isolated contact discontinuity
    for(int i=0 ; i<Nx ; i++){
        if( x[i] < 0.3 ){
            rho[i] = 1.4; // density
        }
        else{
            rho[i] = 1.0; // density
        }

        u[i] = 0.1; // velocity
        p[i] = 1.0; // pressure

        m[i] = rho[i]*u[i]; // momentum
        E[i] = (p[i]/(CTS::GAMMA-1)) + (rho[i]*(u[i]*u[i])/2); // energy
    }

    // 1.2.1  Extend the cells with ghost values on both ends of the conserved variables
    rho.insert(rho.begin(), rho.back());
    rho.push_back(rho[1]);

    m.insert(m.begin(), m.back());
    m.push_back(m[1]);

    E.insert(E.begin(), E.back());
    E.push_back(E[1]);

    // Log
    cout << "Conserved Variables (Initially from Riemann problem) => (len=" << rho.size() << ")" << endl;
    cout << "rho : ";
    for(auto i : rho){
        cout << i << " ";
    }
    cout << endl;
    cout << "m : ";
    for(auto i : m){
        cout << i << " ";
    }
    cout << endl;
    cout << "E : ";
    for(auto i : E){
        cout << i << " ";
    }
    cout << endl;

    cout << "-------------------------------------------------------" << endl;    
    cout << "-------CALCULATING PIECEWISE LINEAR RECONSTRUCTION-----" << endl;
    cout << "-------------------------------------------------------" << endl;

    /*
        STEP-2
        -> Calculation of the Piecewise Linear Reconstruction U_East and U_West
        -> Calculation of Slopes (Ux)j
        -> Extend the cells for the Slope vector
        -> Find the Piecewise Linear Reconstruction
    */
    vector<vector<double>> rho_plr = PLR::construct_plr(rho);
    vector<vector<double>> m_plr = PLR::construct_plr(m);
    vector<vector<double>> E_plr = PLR::construct_plr(E);

    // Log
    PLR::print_plr_conserved(rho_plr, m_plr, E_plr);

    cout << "-------------------------------------------------------" << endl;    
    cout << "---CALCULATING ONE SIDED LOCAL SPEEDS OF PROPAGATION---" << endl;
    cout << "-------------------------------------------------------" << endl;

    /*
        STEP-3 
        -> Calculate the one-sided local speeds of propagation
    */
    // 2D vector to contain the One sided local speeds of propagation a+ and a-
    vector<vector<double>> lsp = LSP::cu_lsp(rho_plr, m_plr, E_plr);
    LSP::print_lsp(lsp); // Log

    //! For the CURH scheme there will be some steps after this 
    //! to modify the values of a+ and a- for better estimates
    //! these modifications are skipped in the CU scheme

    cout << "-------------------------------------------------------" << endl;    
    cout << "------------------CALCULATING FLUX VECTORS-------------" << endl;
    cout << "-------------------------------------------------------" << endl;

    /*
        STEP-4 
        -> Calculate the Anti-diffusion term (from the equation 2.3 in the research paper)
        -> 1) calculate U*
        -> 2) Calculate F(U_West) and F(U_East) for the above purpose
    */

    // 4.1 Calculate the flux vectors
    vector<vector<double>> f1 = FLX::get_flux(1, rho_plr, m_plr, E_plr);
    vector<vector<double>> f2 = FLX::get_flux(2, rho_plr, m_plr, E_plr);
    vector<vector<double>> f3 = FLX::get_flux(3, rho_plr, m_plr, E_plr);

    // log
    cout << "Flux vectors" << endl;

    cout << "f1_E : ";
    for(auto i :f1[0]){
        cout << i << " ";
    }
    cout << endl;
    cout << "f1_W : ";
    for(auto i :f1[1]){
        cout << i << " ";
    }
    cout << endl;

    cout << "f2_E : ";
    for(auto i :f2[0]){
        cout << i << " ";
    }
    cout << endl;
    cout << "f2_W : ";
    for(auto i :f2[1]){
        cout << i << " ";
    }
    cout << endl;

    cout << "f3_E : ";
    for(auto i :f3[0]){
        cout << i << " ";
    }
    cout << endl;
    cout << "f3_W : ";
    for(auto i :f3[1]){
        cout << i << " ";
    }
    cout << endl;

    // 4.2 Calculate U* using a+,a-, U_East, U_West, F_East, F_West
    // 4.3 Calculate dU using U* and other terms
    // 4.4 Calculate CU Numerical Flux using dU and other terms
    vector<vector<double>> star = ADC::CalculateStar(lsp, f1, f2, f3, rho_plr, m_plr, E_plr);
    vector<vector<double>> ADTerm = ADC::CalculateADT(star, rho_plr, m_plr, E_plr);

    vector<vector<double>> cu_flux = CUF::get_cu_flux(lsp, ADTerm, rho_plr, m_plr, E_plr, f1, f2, f3);

    /*
        STEP-5 Calculate the next iteration using SSPRK and CFL
        -> Using CFL Conditions get dT for the next iteration
        -> Now use the 3-stage SSPRK scheme to get the new values for the Conserved Variables
    */

    // 5.1 calculate dT using the CFL conditions
    dT = CFL::get_new_dt(lsp, T_final, dx);

    // 5.2 Using the 3-stage SSPRK scheme we get the new values for the Conserved variables
    


    return 0;
}

