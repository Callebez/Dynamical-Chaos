#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
// #include "../include/gnuplot-iostream.h"
// #include "../include/rungekutta4thSquare.hpp"
#include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"
int main()
{

    // std::vector<double> integrationAux = {6,1,0,1};
    // std::vector<double> auxVec (4,0);
    // double step = 1e-3;
    
    // double timespan[] = {0,10};
    // auxVec = rungeKutta4thSquare(classicalPendulum, integrationAux, 1e-3, step, 4);
    // completeRungeKuttaToFile(classicalPendulum, integrationAux, 1, step, 4, timespan);

    std::vector<std::vector<double>> Matrix;  
    std::vector<double> line (4, 0);
    std::vector<double> line2 (4, 1);
    Matrix.push_back(line);
    Matrix.push_back(line2);

    for(uint i = 0; i < Matrix.size(); i++)
    {
        for(uint j = 0; j < Matrix[0].size(); j++)
        {
            std::cout<<Matrix[i][j] <<" ";
        }
        std::cout<<std::endl;

    }
    return 0;    
}
