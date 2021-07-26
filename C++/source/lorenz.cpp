#include "../include/lorenz.hpp"
/*
    lorenz:     Set of differential equations for the lorenz system, returns a point
    coord:      Point where the system is acting on 
    rho:        Parameter for the biffurcation
*/
std::vector<double> lorenz(std::vector<double> coord, double rho)
{
    std::vector<double> xcoord (3,0);
    
    xcoord[0] = 10.0*coord[1] - 10.0 * coord[0];
    xcoord[1] = rho*coord[0] - coord[1] - coord[0]*coord[2];
    xcoord[2] = coord[0]*coord[1] -(8.0/3.0)*coord[2];
    return xcoord;
}

/*
    lorenzJacobian:     Jacobian for the Lorenz system: J*step + Identity Matrix
    coord:              Point where the system is acting on
    step:               Step of the "line" for which the jacobian is the coefficient 
*/
std::vector<std::vector<double>> lorenzJacobian(std::vector<double>& coord,double step)
{
    std::vector<std::vector<double>> jacobian (3,std::vector<double>(3,0));

    jacobian[0][0] = -10.0*step+1.0; 
    jacobian[0][1] =  10.0*step;
    jacobian[0][2] =  0.0;  
    jacobian[1][0] = -coord[2]*step + 28.0*step;
    jacobian[1][1] = -step+1.0;
    jacobian[1][2] = -coord[0]*step;
    jacobian[2][0] = coord[1]*step;
    jacobian[2][1] = coord[0]*step;
    jacobian[2][2] = 1.0-(8.0/3.0)*step;
    return jacobian;
}