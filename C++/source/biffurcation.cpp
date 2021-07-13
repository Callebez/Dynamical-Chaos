#include "../include/biffurcation.hpp"


std::vector<std::vector<double>> biffurcation(std::vector<double>(*function)(std::vector<double>, double),
                                            double paramRange[2], std::vector<double> initialCond, 
                                            double paramStep, double integrationStep, int systemDimension)
{

    int paramIterations = (int)((fabs(paramRange[1]-paramRange[0])/paramStep));
    int integrationIterations = (int)((double)10.0/integrationStep);
    std::vector<std::vector<double>> biffurcation (4);
    std::vector<std::vector<double>> auxCoord(integrationIterations);
    auxCoord[0]= initialCond;

    std::vector<double> param (paramIterations,0);
    std::vector<double>::iterator paramValue;
    std::vector<double> xCoord (integrationIterations-1,0);
    std::vector<double> auxIntegration;
    for(int i = 0; i < paramIterations; i++)
    {
        param[i] = i*paramStep;
        // std::cout<<param[i]<<" \n";
    } 
    
    
    for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    {
        //integration for transient state
        for(int i = 0; i < integrationIterations-1; i++)
        {
            // std::cout<<auxCoord[i][0];
            auxCoord[i+1] = rungeKutta4thSquare(function, auxCoord[i], *paramValue, integrationStep,systemDimension);
            xCoord[i] = auxCoord[i][2];
           // std::cout<<xCoord[i]<<"\n";
        }
        // for(uint i = 0; i < auxCoord.size(); i++)
        // {
        //     for(uint j = 0; j < auxCoord[0].size(); j++)
        //     {
        //         std::cout<<auxCoord[i][j] <<" ";
        //     }
        //     std::cout<<"\n";
        // }
        //     std::cout<< "   FIM     "<<std::endl;

        auxCoord[0] = auxCoord[integrationIterations-1];
        // std::cout<<auxIntegration[0]<<"\n";
        
     
        //  for(uint i = 0; i < auxCoord.size(); i++)
        // {
        //     for(uint j = 0; j < auxCoord[0].size(); j++)
        //     {
        //         std::cout<<auxCoord[i][j] <<" ";
        //     }
        //     std::cout<<"\n";
        // }
        for(int i = 1; i < integrationIterations-2; i++)
        {
            if( (xCoord[i-1]<xCoord[i]) & (xCoord[i]>xCoord[i+1]))
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
    
