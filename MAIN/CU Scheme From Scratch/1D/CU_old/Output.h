#pragma once

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<utility>

#include "Constants.h"
#include "Parameters.h"

using namespace std;

namespace CTS = Constants;
namespace PRM = Parameters;

#define ll long long int
#define vd vector<double>
#define vvd vector< vector<double> >

namespace Output{

    // Function to convert a vector to a string of space separated values
    string vector_to_string(vd grid);

    // Function to write the computational grid to output file
    void write_computational_grid(vd grid);

    // Function to write the conserved vectors to their output files
    void write_conserved_variables(vvd cons_vars);
}

string Output::vector_to_string(vd grid){
    std::stringstream ss;

    for(auto i : grid){
        ss << i << " ";
    }

    return ss.str();    
}

void Output::write_computational_grid(vd grid){
    string file_name;

    switch (PRM::MODE){
        case 1: // CU
            file_name = "CU-RESULTS/ComputationalDomain.txt";

            break;

        case 2: // CU-REF
            file_name = "CU-REF-RESULTS/ComputationalDomain.txt";

            break;

        case 3: // CURH
            file_name = "CURH-RESULTS/ComputationalDomain.txt";

            break;

        case 4: // CURH-REF
            file_name = "CURH-REF-RESULTS/ComputationalDomain.txt";
             
            break;

        default:
            cout << "---ERROR--- Please enter correct mode ---" << endl;
            break;
    }

    // convert the vectors to strings
    string str_grid = vector_to_string(grid);

    // open required file in append mode
    ofstream fout;
    fout.open(file_name, ios::app);

    // write to file
    fout << str_grid << endl;

    fout.close();

    cout << "---LOG--- write successfull in [" << file_name << "] ---" << endl; 
}

void Output::write_conserved_variables(vvd cons_vars){
    string file_name_1;
    string file_name_2;
    string file_name_3;

    switch (PRM::MODE){
        case 1: // CU

            file_name_1 = "CU-RESULTS/density.txt";
            file_name_2 = "CU-RESULTS/momentum.txt";
            file_name_3 = "CU-RESULTS/energy.txt";

            break;

        case 2: // CU-REF

            file_name_1 = "CU-REF-RESULTS/density.txt";
            file_name_2 = "CU-REF-RESULTS/momentum.txt";
            file_name_3 = "CU-REF-RESULTS/energy.txt";

            break;

        case 3: // CURH

            file_name_1 = "CURH-RESULTS/density.txt";
            file_name_2 = "CURH-RESULTS/momentum.txt";
            file_name_3 = "CURH-RESULTS/energy.txt";

            break;

        case 4: // CURH-REF

            file_name_1 = "CURH-REF-RESULTS/density.txt";
            file_name_2 = "CURH-REF-RESULTS/momentum.txt";
            file_name_3 = "CURH-REF-RESULTS/energy.txt";
            
            break;

        default:
            cout << "---ERROR--- Please enter correct mode ---" << endl;
            break;
    }    

    // convert the vectors to strings
    string str_density = vector_to_string(cons_vars[0]);
    string str_momentum = vector_to_string(cons_vars[1]);
    string str_energy = vector_to_string(cons_vars[2]);

    // open required file in append mode
    ofstream rho_out;
    ofstream m_out;
    ofstream e_out;
    
    rho_out.open(file_name_1, ios::app);
    m_out.open(file_name_2, ios::app);
    e_out.open(file_name_3, ios::app);

    // write to file
    rho_out << str_density << endl;
    m_out << str_momentum << endl;
    e_out << str_energy << endl;

    rho_out.close();
    m_out.close();
    e_out.close();

    cout << "---LOG--- write successfull in [" << file_name_1 << "] ---" << endl; 
    cout << "---LOG--- write successfull in [" << file_name_2 << "] ---" << endl; 
    cout << "---LOG--- write successfull in [" << file_name_3 << "] ---" << endl; 
}
