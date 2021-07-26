#include "../include/biffurcation.hpp"
/* 
    biffurcation:       Function for the evaluation of the biffurcations in the system
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

std::vector<std::vector<double>> biffurcation(std::vector<double>(*function)(std::vector<double>, double),
                                            double paramRange[2], std::vector<double> initialCond, 
                                            double paramStep, double integrationStep, int coordBeingAnalysed  ,int systemDimension)
{

    int paramIterations = (int)((fabs(paramRange[1]-paramRange[0])/paramStep));    
    int integrationIterations = (int)((double)10.0/integrationStep);

    std::vector<std::vector<double>> biffurcation (4);
    std::vector<std::vector<double>> auxCoord(integrationIterations);
    //inicialization for the integration
    auxCoord[0]= initialCond;

    std::vector<double> param (paramIterations,0);
    std::vector<double>::iterator paramValue;
    std::vector<double> xCoord (integrationIterations-1,0);
    std::vector<double> auxIntegration;

    for(int i = 0; i < paramIterations; i++)
    {
        param[i] = paramRange[0] + i*paramStep;
        // std::cout<<param[i]<<" \n";
    } 
    
    
    for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    {
        //integration for transient state
        for(int i = 0; i < integrationIterations-1; i++)
        {
            auxCoord[i+1] = rungeKutta4thSquare(function, auxCoord[i], *paramValue, integrationStep,systemDimension);
            xCoord[i] = auxCoord[i][coordBeingAnalysed];
        }

        auxCoord[0] = auxCoord[integrationIterations-1];

        for(int i = 1; i < integrationIterations-2; i++)
        {
            if((xCoord[i-1]<xCoord[i]) & (xCoord[i]>xCoord[i+1]))
            {
                biffurcation[1].push_back(xCoord[i]);
                biffurcation[0].push_back((double)(*paramValue));
                if(xCoord[i+1]<xCoord[i+2])
                {
                    biffurcation[2].push_back(xCoord[i+1]);
                    biffurcation[3].push_back((double)(*(paramValue + 1)));
                    i++;
                }
            }
            else if((xCoord[i]<xCoord[i-1]) & (xCoord[i]< xCoord[i+1]))
            {
                biffurcation[3].push_back(xCoord[i]);
                biffurcation[2].push_back((double)(*paramValue));
                if(xCoord[i+1]>xCoord[i+2])
                {
                    biffurcation[0].push_back(xCoord[i+1]);
                    biffurcation[1].push_back((double)(*(paramValue + 1)));
                    i++;
                }
            }
        }
        xCoord.clear();
           
      //  std::cout<< *paramValue<<" \n";

    }
    return biffurcation;
}
/*
    printBiffucationToFile:         Prints the matrix, result of the function "biffurcation"
                                    into two files, one containing the maxima and another one 
                                    contaning the minima points as a result of the process of 
                                    biffurcation.
    fileName:                       Name of the file in which the matrix will be printed at

*/
void printBiffucationToFile(std::vector<std::vector<double>>& matrix, std::string fileName)
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
