#include "../include/biffurcation.hpp"


std::vector<std::vector<double>> biffurcation(std::vector<double>(*function)(std::vector<double>, double),
                                            double paramRange[2], std::vector<double> initialCond, 
                                            double paramStep, double integrationStep, int systemDimension)
{
    int paramIterations = (int)(abs(paramRange[1]-paramRange[0])/paramStep);
    int integrationIterations = (int)(20/integrationStep);
    std::vector<std::vector<double>> biffurcation (4);


    std::vector<double> param (paramIterations,0);
    std::vector<double>::iterator paramValue;
    std::vector<double> xCoord (integrationIterations,0);
    for(int i = 0; i < paramIterations; i++)
    {
        param[i] = i*paramStep;
    } 
    
    
    for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    {
        //integration for transient state
        for(int i = 0; i < integrationIterations; i++)
        {
            initialCond = rungeKutta4thSquare(function, initialCond, *paramValue, integrationStep,systemDimension);
            xCoord[i] = initialCond[0];
        }
        for(uint i = 1; i < integrationIterations-1; i++)
        {
            if(xCoord[i-1]<xCoord[i] && xCoord[i]>xCoord[i+1])
            {
                biffurcation[1].push_back(xCoord[i]);
                biffurcation[0].push_back(*paramValue);
                // if(xCoord[i+1]<xCoord[i+2])
                // {
                //     biffurcation[2].push_back(xCoord[i+1]);
                //     biffurcation[3].push_back(*(paramValue + 1));
                //     i++;
                // }
            }
            else if(xCoord[i]<xCoord[i-1] && xCoord[i]< xCoord[i+1])
            {
                biffurcation[3].push_back(xCoord[i]);
                biffurcation[2].push_back(*paramValue);
                // if(xCoord[i+1]>xCoord[i+2])
                // {
                //     biffurcation[0].push_back(xCoord[i+1]);
                //     biffurcation[1].push_back(*(paramValue + 1));
                //     i++;
                // }
            }
        }
    }
        return biffurcation;
}
    
