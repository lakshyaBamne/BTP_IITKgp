/*
    * @author - lakshya Bamne (20MA20029)
    * @supervisor - Prof. Naveen Kumar Garg

    ! Implementation of CU scheme for modified Flux functions in Euler's Equations

    -> Research Paper - 
    Two dimensional Riemann problem for 2x2 system of hyperbolic 
    conservation laws involving three constant states 
    - by Jinah Hwang, Myoungin Shin, Suyeion Shin and Woonjae Hwang
*/
#include<iostream>
#include<string>

#include "EulerGas_2D.h"

using namespace std;

// function to run the Central Upwind scheme for one problem
void RunCentralUpwindScheme(int shock_type, string mode){

    // initialize Riemann problem using the class constructor
    ModifiedEulerSystem_CU rp(shock_type);

    // rp.show_cons_vars_new();

    if( mode == "COMPLETE" ){
        rp.RunCU_complete();
    }
    else if( mode == "PARTIAL" ){
        rp.RunCU_partial();
    }
    else{
        cout << "---ERROR--- Please enter correct mode(COMPLETE/PARTIAL ---" << endl;
    }

}

int main(){
    
    // repeatedly run the Central Upwind Scheme for all initial conditions and store the results
    // for(int i=1 ; i<=12 ; i++){
    //     RunCentralUpwindScheme(i, "PARTIAL");
    // }

    RunCentralUpwindScheme(1, "PARTIAL");

    return  0;
}


