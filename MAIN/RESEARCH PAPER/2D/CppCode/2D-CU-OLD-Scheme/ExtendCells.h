#pragma once

#include<iostream>
#include<vector>
#include<utility>

using namespace std;

namespace ExtendCells{
    //! 1D code
    void extend_cells(string bound_cond, vector<vector<double>>& cons_vars);
    void extend_cells_ind(string bound_cond, vector<double>& cons_vars, string var);
    void extend_plr(string bound_cond, string var, vector<double>& east, vector<double>& west);

    //! 2D code
    void extend_cells_horizontal(string bound_cond, vector< vector<double> >& cons_vars);
    void extend_cells_vertical(string bound_cond, vector< vector<double> >& cons_vars);
    void extend_plr_horizontal(string bound_cond, string var, vector<double>& east, vector<double>& west);
    void extend_plr_vertical(string bound_cond, string var, vector<double>& east, vector<double>& west);
}

//! 2D
void ExtendCells::extend_cells_horizontal(string bound_cond, vector< vector<double> >& cons_vars){
    if( bound_cond == "FREE" ){

        for(int u=0 ; u<4 ; u++){
            cons_vars[u].insert(cons_vars[u].begin(), cons_vars[u][0]);
            cons_vars[u].push_back(cons_vars[u].back());
        }

    }
    else if( bound_cond == "RFREE" ){

        for(int u=0 ; u<4 ; u++){
        
            if( u!=1 ){
                cons_vars[u].insert(cons_vars[u].begin(), cons_vars[u][0]);
                cons_vars[u].push_back(cons_vars[u].back());
            }
            else{
                // extend m1
                cons_vars[u].insert(cons_vars[u].begin(), -cons_vars[u][0]);
                cons_vars[u].push_back(cons_vars[u].back());
            }
        
        }

    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}

void ExtendCells::extend_cells_vertical(string bound_cond, vector< vector<double> >& cons_vars){
    if( bound_cond == "FREE" ){

        for(int u=0 ; u<4 ; u++){
            cons_vars[u].insert(cons_vars[u].begin(), cons_vars[u][0]);
            cons_vars[u].push_back(cons_vars[u].back());
        }

    }
    else if( bound_cond == "RFREE" ){

        for(int u=0 ; u<4 ; u++){
        
            if( u!=2 ){
                cons_vars[u].insert(cons_vars[u].begin(), cons_vars[u][0]);
                cons_vars[u].push_back(cons_vars[u].back());
            }
            else{
                // extend m1
                cons_vars[u].insert(cons_vars[u].begin(), -cons_vars[u][0]);
                cons_vars[u].push_back(cons_vars[u].back());
            }
        
        }

    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}

void ExtendCells::extend_plr_horizontal(string bound_cond, string var, vector<double>& east, vector<double>& west){
    if( bound_cond == "FREE" ){
        east.insert(east.begin(), west[0]);
        west.push_back(east.back());

        east.push_back(0); //! unnecessary
        west.insert(west.begin(),0); //! unnecessary
    }
    else if( bound_cond == "RFREE" ){
        if( var=="momentumX" ){
            east.insert(east.begin(), -west[0]);
            west.push_back(east.back());

            east.push_back(0); //! unnecessary
            west.insert(west.begin(), 0); //! unnecessary            
        }
        else{
            east.insert(east.begin(), west[0]);
            west.push_back(east.back());

            east.push_back(0); //! unnecessary
            west.insert(west.begin(), 0); //! unnecessary
        }
    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}

void ExtendCells::extend_plr_vertical(string bound_cond, string var, vector<double>& east, vector<double>& west){
    if( bound_cond == "FREE" ){
        east.insert(east.begin(), west[0]);
        west.push_back(east.back());

        east.push_back(0); //! unnecessary
        west.insert(west.begin(),0); //! unnecessary
    }
    else if( bound_cond == "RFREE" ){
        if( var=="momentumY" ){
            east.insert(east.begin(), -west[0]);
            west.push_back(east.back());

            east.push_back(0); //! unnecessary
            west.insert(west.begin(), 0); //! unnecessary            
        }
        else{
            east.insert(east.begin(), west[0]);
            west.push_back(east.back());

            east.push_back(0); //! unnecessary
            west.insert(west.begin(), 0); //! unnecessary
        }
    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}


//! 1D
void ExtendCells::extend_cells(string bound_cond, vector<vector<double>>& cons_vars){
    if( bound_cond == "FREE" ){

        for(int u=0 ; u<3 ; u++){
            cons_vars[u].insert(cons_vars[u].begin(), cons_vars[u][0]);
            cons_vars[u].push_back(cons_vars[u].back());
        }

    }
    else if( bound_cond == "REFLECTIVE" ){

        for(int u=0 ; u<3 ; u++){
        
            if( u!=1 ){
                cons_vars[u].insert(cons_vars[u].begin(), cons_vars[u][0]);
                cons_vars[u].push_back(cons_vars[u].back());
            }
            else{
                // extend m
                cons_vars[1].insert(cons_vars[1].begin(), -cons_vars[1][0]);
                cons_vars[1].push_back(-cons_vars[1].back());
            }
        
        }

    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}

void ExtendCells::extend_cells_ind(string bound_cond, vector<double>& cons_vars, string var){
    if( bound_cond == "FREE" ){
        cons_vars.insert(cons_vars.begin(), cons_vars[0]);
        cons_vars.push_back(cons_vars.back());
    }
    else if( bound_cond == "REFLECTIVE" ){
        if( var=="momentum" ){
            cons_vars.insert(cons_vars.begin(), -cons_vars[0]);
            cons_vars.push_back(-cons_vars.back());
        }
        else{
            cons_vars.insert(cons_vars.begin(), cons_vars[0]);
            cons_vars.push_back(cons_vars.back());
        }
    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}

void ExtendCells::extend_plr(string bound_cond, string var, vector<double>& east, vector<double>& west){
    if( bound_cond == "FREE" ){
        east.insert(east.begin(), west[0]);
        west.push_back(east.back());

        east.push_back(0); //! unnecessary
        west.insert(west.begin(),0); //! unnecessary
    }
    else if( bound_cond == "REFLECTIVE" ){
        if( var=="momentum" ){
            east.insert(east.begin(), -west[0]);
            west.push_back(-east.back());

            east.push_back(0); //! unnecessary
            west.insert(west.begin(), 0); //! unnecessary            
        }
        else{
            east.insert(east.begin(), west[0]);
            west.push_back(east.back());

            east.push_back(0); //! unnecessary
            west.insert(west.begin(), 0); //! unnecessary
        }
    }
    else{
        cout << "---ERROR--- Please enter correct Boundary Conditions ---" << endl;
    }
}



