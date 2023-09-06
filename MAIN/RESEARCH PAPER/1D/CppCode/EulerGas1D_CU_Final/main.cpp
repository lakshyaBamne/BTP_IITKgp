#include<iostream>
#include<vector>

#include "1D_CentralUpwind.h"
#include "Utility.h"

namespace UTL = Utility;

using namespace std;

int main(){
    string mode;
    
    cout << "+------------------------------------------------------------------------------------+" << endl;
    cout << "|                                Central Upwind Scheme                               |" << endl;
    cout << "+------------------------------------------------------------------------------------+" << endl;
    cout << "|                                     Enter MODE                                     |" << endl;
    cout << "+------------------------------------------------------------------------------------+" << endl;

    cout << "-> [COMPLETE] Single instance , store data for each time step" << endl;
    cout << "-> [PARTIAL] Single instance , store data only for initial and final steps" << endl;
    cout << "-> [PLOT-PARTIAL] Two instances , store data only for initial and final steps" << endl;
    cout << "-> [ANIMATE-ALL] Two instances , store data for each time step" << endl;

    cout << "+------------------------------------------------------------------------------------+" << endl;

    cout << "-> ";
    cin >> mode;

    if( mode != "" ){
        UTL::export_string(mode); // export mode to a file for later use in plotting
    }
    else{
        cout << "---ERROR--- Enter a correct mode to start iterations ---" << endl;
        return 0;
    }


    if( mode == "COMPLETE" ){
        GetInput I;
        I.run_cu_scheme_complete("result");
    }
    else if( mode == "PARTIAL" ){
        GetInput I;
        I.run_cu_scheme_partial("result");
    }
    else if( mode == "PLOT-PARTIAL" ){
        GetInput I1;
        I1.run_cu_scheme_partial("result1");
    
        GetInput I2;
        I2.run_cu_scheme_partial("result2");
    }
    else if( mode == "ANIMATE-ALL" ){
        GetInput I1;
        I1.run_cu_scheme_complete("result1");

        GetInput I2;
        I2.run_cu_scheme_complete("result2");
    }
    else{
        cout << "---ERROR--- Please enter correct MODE to run simulations ---" << endl;
    }

    return 0;
}