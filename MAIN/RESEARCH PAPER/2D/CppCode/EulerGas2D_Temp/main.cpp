/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#include<iostream>
#include<vector>
#include<string>
#include<time.h>
#include<iomanip>

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

    time_t start, end;

    time(&start);

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

    time(&end);

    double time_taken = double(end - start);
    double hours = time_taken/86400;

    cout << "----------------------------------------------------------------------------------------------------" << endl;
    cout << "Time taken for execution : " << fixed << time_taken << setprecision(5);
    cout << "sec" << endl;
    cout << "----------------------------------------------------------------------------------------------------" << endl;
    cout << "Time taken for execution : " << fixed << hours << setprecision(5);
    cout << "hours" << endl;
    cout << "----------------------------------------------------------------------------------------------------" << endl;

    return 0;
}

