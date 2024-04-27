#pragma once

#include "Constants.h"

namespace CTS = Constants;

namespace OutputResult{
    // function to write the computational grid to external file
    void write_grid(const vd& grid, string path);

    // function to write the density to external file
    void write_density(const vd& grid, string path);
}

void OutputResult::write_grid(const vd& grid, string path){
    stringstream ss;
    for(int i=1 ; i<grid.size()-1 ; i++){
        ss << grid[i] << " ";
    }
    
    ofstream fout;
    fout.open(path, ios::app);

    fout << ss.str() << endl;

    fout.close();
}

void OutputResult::write_density(const vd& density, string path){
    stringstream ss;
    for(int i=1 ; i<density.size()-1 ; i++){
        ss << density[i] << " ";
    }

    ofstream fout;
    fout.open(path, ios::app);

    fout << ss.str() << endl;

    fout.close();
}


