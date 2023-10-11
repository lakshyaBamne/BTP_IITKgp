/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#include<iostream>
#include<vector>
#include<string>

#include "Constants.h"
#include "Utility.h"
#include "2D_CUScheme.h"

using namespace std;

namespace CTS = Constants;
namespace CTR = ConstantsRT;
namespace UTL = Utility;

int main(){
    string to_solve;
    string plot_ref;

    cout << "----------------------------------------------------------------------------------------------------" << endl;
    cout << "--------------------------------------------2D Central Upwind---------------------------------------" << endl;
    cout << "----------------------------------------------------------------------------------------------------" << endl;
    
    cout << "---------------------------------------AVAILABLE PROBLEMS---------------------------------------------" << endl;
    cout << "[MCW] Moving Contact Wave" << endl;
    cout << "[2DR] 2D Riemann Problem" << endl;
    cout << "[EXP] Explosion problem" << endl;
    cout << "[IMP] Implosion problem" << endl;
    cout << "[KHI] Kelvin Helmholtz Instability" << endl;
    cout << "[MCW] Rayleigh Taylor Instability" << endl;

    cout << "----------------------------------------------------------------------------------------------------" << endl;
    cout << "SELECT PROBLEM -> ";
    cin >> to_solve;
    cout << "PLOT REFERENCE? (Y/N) -> ";
    cin >> plot_ref;

    if( plot_ref == "Y" ){ // Plot CU and REF both
        EulerSystem_2D P1("CU", to_solve);
        EulerSystem_2D P2("REF", to_solve);

        P1.RunCU_partial();
        P2.RunCU_partial();
    }
    else{ // Plot only CU
        EulerSystem_2D P1("CU", to_solve);

        P1.RunCU_partial();
    }

    return 0;
}

