#include<iostream>

#include "Constants.h"
#include "Initialize.h"
#include "OutputResult.h"
#include "PolynomialReconstruction.h"
#include "InterfaceVariables.h"
#include "LLF_flux.h"
#include "TimeIntegration.h"
#include "TimeStep.h"

using namespace std;

namespace CTS = Constants;
namespace IRP = Initialize;
namespace OPR = OutputResult;
namespace SOR = SecondOrderReconstructions;
namespace IFV = InterfaceVariables;
namespace LLF = LLF;
namespace TIS = TimeIntegrationSchemes;
namespace TUS = TimeStepCalculation;

void update_conserved_variables(vvd& U, vvd& P, vd& a, double& dt);

int main(){

    string MODE = "NORMAL";
    // string MODE = "REF";

    //! STEP-0 : Initialize the variables on the computational grid
    //! --------------------------------------------------------------------------------
    // Initialize the computational grid
    vd grid = IRP::make_grid();

    //! write computational grid to external file
    OPR::write_grid(grid, MODE+"grid.txt");

    // Initialize the primitive variables according to the problem
    vvd P = IRP::set_primitive_variables(grid);

    // Initialize the conserved variables along with other variables
    vvd U(3, vd(CTS::N+2));
    vd a(CTS::N+2);

    IRP::initialize_conserved_variables(P, U, a);

    double dt;
    double t=CTS::time.first;

    while(t < CTS::time.second){
        // log time
        cout << "t=" << t << " | dt=" << dt << endl;

        //! write output to external file
        OPR::write_density(U[0], MODE+"density.txt");

        // update the conserved variables
        update_conserved_variables(U, P, a, dt);

        // update time
        t += dt;
    }

}

void update_conserved_variables(vvd& U, vvd& P, vd& a, double& dt){
    
    //! STEP-1 : Calculate time step
    //! --------------------------------------------------------------------------------
    TUS::update_time(dt, P[1], a);

    //! STEP-2 : Calculate the 2nd order reconstructions
    //! --------------------------------------------------------------------------------
    vvd slx = SOR::get_slopes(U);

    vvvd recons = SOR::get_reconstructions(U, slx);
    vvd Up = recons[0];
    vvd Um = recons[1];

    //! STEP-3 : Calculate Other variables at the interfaces using Up and Um
    //! --------------------------------------------------------------------------------
    vvvd interface_fluxes = IFV::get_interface_fluxes(Up, Um);
    vvd Fp = interface_fluxes[0];
    vvd Fm = interface_fluxes[1];

    vvvd prim_vars = IFV::get_primitive_variables(Up, Um);
    vvd Pp = prim_vars[0];
    vvd Pm = prim_vars[1];

    vvd acos = IFV::get_acoustic_speeds(Pp, Pm);

    //! STEP-4 : Calculate intermediate alpha term for the LLF Flux
    //! --------------------------------------------------------------------------------
    vd alpha = LLF::get_alpha(Pp, Pm, acos);

    //! STEP-5 : Calculate LLF Flux
    //! --------------------------------------------------------------------------------
    vvd llf_flx = LLF::get_llf_fluxes(Fp, Fm, alpha, Up, Um);

    //! STEP-6 : Use time discretization for calculating U(n+1)
    //! --------------------------------------------------------------------------------
    TIS::EFD(U, llf_flx, dt);
    TIS::update_variables(U, P, a);
}
