#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
// #include "../include/gnuplot-iostream.h"
 //#include "../include/rungekutta4thSquare.hpp"
//#include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"
#include "../include/LyapExp.hpp"
#include "../include/lorenz.hpp"

int main()
{
    std::vector<double> integrationAux = {1.0,2.0,3.0};

    std::vector<long double> lya = lyapunovSpectrum(lorenz,lorenzJacobian,integrationAux,0.001,28.0);
    // std::vector<std::vector<double>> A = matMult(M,N);
    std::cout<<"lyapunov exponents: "<<(lya[0])<<", "<<(lya[1])<<", "<<(lya[2])<<"\n ";
    std::cout<<"lyapunov numbers: "<<exp(lya[0])<<", "<<exp(lya[1])<<", "<<exp(lya[2])<<"\n ";

    std::cout<<"sum of lyapunov exponents: "<< lya[0] + lya[1] + lya[2];
    return 0;    
}
