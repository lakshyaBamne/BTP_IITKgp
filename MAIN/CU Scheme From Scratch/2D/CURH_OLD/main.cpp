/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#include<iostream>
#include<vector>
#include<string>

#include "2D_CUScheme.h"

using namespace std;

int main(){

    string to_solve;
    string only_cu;

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
    cin >> only_cu;

    // class RP (RiemannProblem) is responsible for initializing variables and running Numerical Schemes
    if( only_cu == "Y" ){
        // plot both the CU and Reference plot
        RP rp1("CU", to_solve);
        RP rp2("REF", to_solve);

        rp1.RunCU_partial();
        rp2.RunCU_partial();
    }
    else{
        // plot only the CU plot
        RP rp1("CU", to_solve);

        rp1.RunCU_partial();
    }

    return 0;
}

