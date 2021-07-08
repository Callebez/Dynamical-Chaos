#include "../include/biffurcation.hpp"


std::vector<std::vector<double>> biffurcation(std::vector<double>(*function)(std::vector<double>, double),
                                            double paramRange[2], std::vector<double> initialCond, 
                                            double paramStep, double integrationStep, int dimension)
{
    int paramIterations = (int)(abs(paramRange[1]-paramRange[0])/paramStep);
    std::vector<std::vector<double>> coords (2000+1);
    std::vector<std::vector<double>> results;
    std::vector<double> param (paramIterations,0);
    std::vector<double>::iterator ptr;
    std::vector<double> x_maxes;
    std::vector<double> x_mins;
    std::vector<double> param_maxes;
    std::vector<double> param_mins;
    for ( int i = 0 ; i < 2000+1 ; i++ ){coords[i].resize(dimension);}

   
    for(int i = 0; i < paramIterations; i++)
    {
        param[i] = i*paramStep;
    } 
    for(ptr = param.begin(); ptr < param.end(); ptr++)
    {
        for(int i = 0; i < 2000; i++) // Transiente evolution 
        {
            coords[i+1] = rungeKutta4thSquare(function,coords[i], *ptr, integrationStep, dimension);
        }
        for(int i = 1; i < 2000-1; i++)
        {
            if(coords[i-1][0] < coords[i][0]  && coords[i][0] > coords[i + 1][0])
            {
                x_maxes.push_back(coords[i][0]);
                param_maxes.push_back(*ptr);
            }
            if(coords[i-1][0] > coords[i][0]  && coords[i][0] < coords[i + 1][0])
            {
                x_mins.push_back(coords[i][0]);
                param_mins.push_back(*ptr);
            }
        }
        coords[0] = coords[2001];
    }
    results.push_back(x_maxes);
    results.push_back(x_mins);
    results.push_back(param_maxes);
    results.push_back(param_mins);
    return results;
}