/*
    ! Central Upwind Scheme (Old version)
    @author : Lakshya Bamne (20MA20029), 4th Year UG at IIT Kharagpur
    @supervisor : Dr. Naveen Kumar Garg (Department of Mathematics)
*/

#include<iostream>
#include<vector>

#include "1D_CentralUpwind.h"
#include "Utility.h"

namespace UTL = Utility;

using namespace std;

int main(){
    string mode;

    cout << "Central Upwind Scheme (old scheme)" << endl;
    cout << "Enter MODE" << endl;

    cout << "-> [PARTIAL] Single instance , store data only for initial and final steps" << endl;
    cout << "-> [PLOT-PARTIAL] Two instances , store data only for initial and final steps" << endl;

    cout << "-> ";
    cin >> mode;

    if( mode != "" ){
        UTL::export_string(mode); // export mode to a file for later use in plotting
    }
    else{
        cout << "---ERROR--- Enter a correct mode to start iterations ---" << endl;
        return 0;
    }

    if( mode == "PARTIAL" ){
        GetInput I;
        I.run_cu_scheme_partial("result");
    }
    else if( mode == "PLOT-PARTIAL" ){ // for plotting result graphs with reference
        GetInput I1;
        I1.run_cu_scheme_partial("result1");
    
        GetInput I2;
        I2.run_cu_scheme_partial("result2");
    }
    else{
        cout << "---ERROR--- Please enter correct MODE to run simulations ---" << endl;
    }

    return 0;
}