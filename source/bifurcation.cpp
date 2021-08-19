#include "../include/bifurcation.hpp"
#include "../include/printing.hpp"
#include <thread>

/* 
    bifurcation:       Function for the evaluation of the biffurcations in the system
                        It works by varing the parameter and integrating the system of 
                        equations each time, and them searching the coordinate for points
                        where the function assumes a maxima ou a minima

    function:           Set of differential equations
    paramRange          Interval in which the select parameter will vary
    initalCond          Initial Conditions for the integrations of the system
    paramStep           Step size in between each parameter integration
    integrationStep     Step size for the rk4 integration
    coordBeingAnalysed  Selects in which coordinate the biffucations are being analysed
    systemDimension     Dimension of the system of equations
*/

std::vector<std::vector<double>> bifurcation(std::vector<double>(*function)(std::vector<double>, double),
                                            double paramRange[2], std::vector<double> initialCond, 
                                            double paramStep, int coordBeingAnalysed, int systemDimension)
{

    int paramIterations = (int)((fabs(paramRange[1]-paramRange[0])/paramStep));    

    std::vector<std::vector<double>> biffurcation (4);
   //inicialization for the integration

    std::vector<double> param (paramIterations,0);
    std::vector<double>::iterator paramValue;
    std::vector<double> xCoord;
    double integrationStep = 0.1;
    for(int i = 0; i < paramIterations; i++)
    {
        param[i] = paramRange[0] + i*paramStep;
    } 
    std::vector<std::vector<double>> rk45;
    std::vector<double> hs;
   

    for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    {
        completeRungeKutta45(function,initialCond, rk45, *paramValue,integrationStep,initialCond.size(),1e-3, 100, hs);

        xCoord.resize(rk45.size());
        for(uint i = 0; i< rk45.size(); i++)
        {
            xCoord[i] = rk45[i][coordBeingAnalysed];
        }
 
        initialCond = rk45.back();
        rk45.resize(0);
 
        searchMaxMin(biffurcation,xCoord, paramValue);           

    }
    return biffurcation;
}
void searchMaxMin(std::vector<std::vector<double>>& biffurcation,std::vector<double>& xCoord, std::vector<double>::iterator& paramValue)
{

        for(int i = 0; i < (int)xCoord.size()-1; i++)
        {
            if((xCoord[i-1]<xCoord[i]) & (xCoord[i]>xCoord[i+1]))
            {
                biffurcation[1].push_back(xCoord[i]);
                biffurcation[0].push_back((double)(*paramValue));
            }
            else if((xCoord[i]<xCoord[i-1]) & (xCoord[i]< xCoord[i+1]))
            {
                biffurcation[3].push_back(xCoord[i]);
                biffurcation[2].push_back((double)(*paramValue));
            }
        }
}
/*
    printBifucationToFile:         Prints the matrix, result of the function "bifurcation"
                                    into two files, one containing the maxima and another one 
                                    contaning the minima points as a result of the process of 
                                    biffurcation.
    fileName:                       Name of the file in which the matrix will be printed at

*/
void printBifucationToFile(std::vector<std::vector<double>>& matrix, std::string fileName)
{
    std::ofstream output1;
    std::ofstream output2;

    std::string file1 = "./outputs/txt/";
    file1 = file1 + fileName + "max.dat";
    
    std::string file2 = "./outputs/txt/";
    file2 = file2 + fileName + "min.dat";
    output1.open(file1);
    output2.open(file2);
    for(uint i = 0; i < 2; i++)
    {
        for(uint j = 0; j < matrix[0].size(); j++)
        {
            output1<<matrix[0][j]<<" "<< matrix[1][j] <<"\n";
            output2<<matrix[2][j]<<" "<< matrix[3][j] <<"\n";

        }
        output1<<" ";
        output2<<" ";

    }
    output1.close();
    output2.close();
}
