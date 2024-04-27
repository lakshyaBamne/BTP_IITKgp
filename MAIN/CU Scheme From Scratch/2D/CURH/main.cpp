/*
    ! @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-8 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme - with RH optimization for Euler Equations of Gas Dynamics (CURH Scheme)
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<unordered_map>
#include<cmath>

#include "Constants.h"
#include "Utility.h"
#include "CURH.h"

namespace CTS=Constants;
namespace UTL=Utility;

using namespace std;

int main(){
    
    //! Create objects for all the required simulations
    // CURH sim1(1,1); //! MCW
    CURH sim2(2,1); //! 2DR-CFG3
    // CURH sim3(3,1); //! EXP
    // CURH sim4(4,1); //! IMP1
    // CURH sim5(5,1); //! IMP2
    // CURH sim6(6,1); //! KHI1
    // CURH sim7(7,1); //! KHI2
    // CURH sim8(8,1); //! KHI3
    // CURH sim9(9,1); //! RTI

    //! Reference simulations with finer computational domain (4x finer)
    // CURH ref1(1,2); //! MCW
    // CURH ref2(2,2); //! 2DR-CFG3
    // CURH ref3(3,2); //! EXP
    // CURH ref4(4,2); //! IMP1
    // CURH ref5(5,2); //! IMP2
    // CURH ref6(6,2); //! KHI1
    // CURH ref7(7,2); //! KHI2
    // CURH ref8(8,2); //! KHI3
    // CURH ref9(9,2); //! RTI

    // sim1.RunCURH_partial();
    sim2.RunCURH_partial();
    // sim3.RunCURH_partial();

    return 0;
}
