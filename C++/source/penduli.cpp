#include "../include/penduli.hpp"
#include <iostream>
std::vector<double> classicalPendulum(std::vector<double> coord, double k)
{
    std::vector<double> coord_dot (4);
    coord_dot[0] = coord[2];
    coord_dot[1] = coord[3];
    coord_dot[2] = -coord[0] - k*coord[0]*coord[1]*coord[1];
    coord_dot[3] = -coord[1] - k*coord[1]*coord[0]*coord[0];
    
    return coord_dot;
}
std::vector<std::vector<double>> classicalPendulumJacobian(std::vector<double>& coord, double k, double step)
{
    std::vector<std::vector<double>> jacobian (4,std::vector<double>(4,0));
    jacobian[0][0] = 1.0;
    jacobian[0][1] = 0.0;
    jacobian[0][2] = step;
    jacobian[0][3] = 0.0;
    
    jacobian[1][0] = 0.0;
    jacobian[1][1] = 1.0;
    jacobian[1][2] = 0.0;
    jacobian[1][3] = step;

    //jacobian[2][0] = -step;
    //jacobian[2][1] = 0.0;
    //jacobian[2][2] = 1.0;
    //jacobian[2][3] = 0.0;
    //
    //jacobian[3][0] = 0.0;
    //jacobian[3][1] = 0.0;
    //jacobian[3][2] = -step;
    //jacobian[3][3] = 1.0;

    jacobian[2][0] = -step-coord[1]*coord[1]*step;
    jacobian[2][1] = -2*coord[0]*coord[1]*step;
    jacobian[2][2] = 1.0;
    jacobian[2][3] = 0.0;
    jacobian[3][0] = -2 * coord[0] * coord[1] * step;
    jacobian[3][1] = -step-coord[0]*coord[0]*step;
    jacobian[3][2] = 0.0;
    jacobian[3][3] = 1.0;

    return jacobian;
}
std::vector<double> quantumPendulum(std::vector<double> coord, double gamma)
{
    std::vector<double> coord_dot (4,0);
    double f = flutuation(coord[0],coord[1], gamma);
    coord_dot[0] = coord[2];
    coord_dot[1] = coord[3];
    coord_dot[2] = - gamma*coord[2] -coord[0] - coord[0]*coord[1]*coord[1]  + f;
    coord_dot[3] = -gamma*coord[3] - coord[1]- coord[0]*coord[0]*coord[1] + f;
 
    return coord_dot; 
}
void waveFunction(double x, double y, double& rho)
{
    double psi = exp(-0.9122350052*(x*x + y*y))*(1 + 0.3345167463*(x*x+y*y)+0.09871707060*(pow(x,4)+pow(y,4)));
    rho = pow(psi,2);
}

double flutuation(double x, double y, double gamma)
{
    srand (time(NULL));
    double a = (float)rand()/RAND_MAX;
    return  a*sqrt((gamma)*((0.62832912000*(pow(x,2.0) + pow(y,2.0)) +
                        0.01205153100* pow(x,2.0) * pow(y,2.0)  +
                        0.56437351570*(pow(x,4.0) + pow(y,4.0)) +
                        0.13998728990*(pow(x,6.0) + pow(y,6.0)) +
                        0.06202680930*(pow(x,2.0)* pow(y,4.0) +  pow(y,2.0)*pow(x,4.0)) +
                        0.03555913955* pow(x,4.0) *  pow(y,4.0)  +
                        0.01777956970*( pow(x,8.0) +  pow(y,8.0)) + 1.155436517)/
                        pow((1.0 + 0.3345167463*(pow(x,2.0) +  pow(y,2.0)) +
                        0.09871707060*(pow(x,4.0) +  pow(y,4.0))),2.0)));
}





