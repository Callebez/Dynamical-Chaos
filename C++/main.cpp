#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
#include "gnuplot-iostream.h"
#include "biffurcation.cpp"
#include "penduli.cpp"

int main()
{

    // std::vector<std::vector<double>> plotPendulum (1000); 
    // for ( int i = 0 ; i < 1000 ; i++ ){plotPendulum[i].resize(4);}
    std::vector<double> integrationAux = {6,1,0,1};
    std::vector<double> auxVec (4,0);
    double step = 1e-4;
    std::ofstream a;
    a.open("teste.dat");
    auxVec = rungeKutta4thSquare(classicalPendulum, integrationAux, 1e-3, step, 4);

    std::ofstream biffplot;
    biffplot.open ("biff.dat");
    for(int i = 0; i < 100000; i++)
    {
        auxVec = rungeKutta4thSquare(classicalPendulum, auxVec, 1e-3, step, 4);
        for(int j = 0; j < 3; j++)
        {
            biffplot << auxVec[j] <<"  ";
        }
        biffplot<< std::endl;
    }
   
 
    // std::vector<double> initalCond = {6,1,0,1};

    // // for(int i = 0; i < 4; i++)
    // // {
    // //     std::cout<< classicalPendulum(initalCond,1)[i]<< "\n";

    // // }
    // // Gnuplot gp;
    // double time[2] = {0,2.0};
    // std::vector<std::vector<double>> biff = biffurcation(classicalPendulum, time, initalCond, 0.1, 0.1,4);
	
    // std::vector<std::vector<double>> print0;
    // print0[0] = biff[0];
    // // print0[1] = biff[1];
    // for(int i = 0; i < 1000; i++)
    // {
    //     for(int j = 0; j < 4; j++)
    //     {
    //         biffplot<< plotPendulum[j][i] << "  ";
    //     }
    //     biffplot<<"\n";
    // }
    // // std::vector<double> x_mins     = biff[1];
    // // std::vector<double> param_maxes= biff[2];
    // // std::vector<double> param_mins = biff[3];

	
	// Data will be sent via a temporary file.  These are erased when you call
	// gp.clearTmpfiles() or when gp goes out of scope.  If you pass a filename
	// (e.g. "gp.file1d(pts, 'mydata.dat')"), then the named file will be created
	// and won't be deleted (this is useful when creating a script).
	// gp << "plot" << gp.file1d(print0) << std::endl;
     biffplot.close();
    a.close();
    return 0;    
}
