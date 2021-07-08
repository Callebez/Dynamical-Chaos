#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
//#include "../include/gnuplot-iostream.h"
#include "../include/rungekutta4thSquare.hpp"
// #include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"
// std::vector<double> classicalPendulum(std::vector<double> coord, double k)
// {
//     std::vector<double> coord_dot (4,0);
//     coord_dot[0] = coord[2];
//     coord_dot[1] = coord[3];
//     coord_dot[2] = -coord[0] - k*coord[0]*coord[1]*coord[1];
//     coord_dot[3] = -coord[1] - k*coord[1]*coord[0]*coord[0];
    
//     return coord_dot;
// }
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
