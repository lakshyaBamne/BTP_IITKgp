/*
    Extend the cells of a given vector of Conserved variables based
    on the given boundary conditions
*/
#pragma once

#include<iostream>
#include<vector>
#include<utility>

using namespace std;

namespace ExtendCells{
    void extend_cells(string bound_cond, vector<vector<double>>& cons_vars);

    void extend_cells_ind(string bound_cond, vector<double>& cons_vars);
}

void ExtendCells::extend_cells(string bound_cond, vector<vector<double>>& cons_vars){
    if( bound_cond == "FREE" ){
        // extend rho
        cons_vars[0].insert(cons_vars[0].begin(), cons_vars[0][0]);
        cons_vars[0].push_back(cons_vars[0].back());

        // extend m
        cons_vars[1].insert(cons_vars[1].begin(), cons_vars[1][0]);
        cons_vars[1].push_back(cons_vars[1].back());

        // extend E
        cons_vars[2].insert(cons_vars[2].begin(), cons_vars[2][0]);
        cons_vars[2].push_back(cons_vars[2].back());
    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}

void ExtendCells::extend_cells_ind(string bound_cond, vector<double>& cons_vars){
    if( bound_cond == "FREE" ){
        cons_vars.insert(cons_vars.begin(), cons_vars[0]);
        cons_vars.push_back(cons_vars.back());
    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}
