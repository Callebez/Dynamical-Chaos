#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
#include "../include/gnuplot-iostream.h"
#include "../include/rungekutta4thSquare.hpp"
#include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"

int main()
{

    std::vector<double> integrationAux = {6,1,0,1};
    std::vector<double> auxVec (4,0);
    double step = 1e-3;
 
    double timespan[] = {0,10};
    auxVec = rungeKutta4thSquare(classicalPendulum, integrationAux, 1e-3, step, 4);
    completeRungeKuttaToFile(classicalPendulum, integrationAux, 1, step, 4, timespan);

    return 0;    
}
