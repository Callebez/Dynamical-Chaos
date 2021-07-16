#include "../include/systems.hpp"

std::vector<double> classicalPendulum(std::vector<double> coord, double k)
{
    std::vector<double> coord_dot (4);
    coord_dot[0] = coord[2];
    coord_dot[1] = coord[3];
    coord_dot[2] = -coord[0] - k*coord[0]*coord[1]*coord[1];
    coord_dot[3] = -coord[1] - k*coord[1]*coord[0]*coord[0];
    
    return coord_dot;
}
std::vector<double> quantumPenduli(std::vector<double> coord, double gamma)
{
    std::vector<double> coord_dot (4,0);
    srand (time(NULL));
    double f = flutuation(coord[0],coord[1], gamma);
    coord_dot[0] = coord[2];
    coord_dot[1] = coord[3];
    coord_dot[2] = - gamma*coord[2] -coord[0] - coord[0]*coord[1]*coord[1] + f*rand();
    coord_dot[3] = -gamma*coord[3] - coord[1]- coord[1]*coord[0]*coord[0]+f*rand();

    return coord_dot; 
}
std::vector<double> lorenz(std::vector<double> coord, double rho)
{
    std::vector<double> coord_dot(3,0);
    coord_dot[0] = 10.0*(coord[1]- coord[0]);
    coord_dot[1] = coord[0]*rho- coord[0]*coord[2] - coord[1];
    coord_dot[2] = coord[0]*coord[1] - 8.0/3.0*coord[2];
    return coord_dot;
    
}

std::vector<std::vector<double>> LorenzJacobian(std::vector<double> coord, double rho)
{
    std::vector<std::vector<double>> jacobian (3);
    double a=10;
    double b=8/3;
    jacobian[0]={-a,a,0};
    jacobian[1]={rho-coord[2],-1,-coord[0]};
    jacobian[2]={coord[1],coord[0],-b};
    return jacobian;
}

double flutuation(double x, double y, double gamma)
{
        return sqrt(
            gamma*(0.62832912000*(pow(x,2.0) + pow(y,2.0)) +
             0.01205153100* pow(x,2.0) * pow(y,2.0)  +
             0.56437351570*(pow(x,4.0) + pow(y,4.0)) +
             0.13998728990*(pow(x,6.0) + pow(y,6.0)) +
             0.06202680930*(pow(x,2.0)* pow(y,4.0) +  pow(y,2.0)*pow(x,4.0)) +
             0.03555913955* pow(x,4.0) *  pow(y,4.0)  +
             0.01777956970*( pow(x,8.0) +  pow(y,8.0)) + 1.155436517)/
             pow((1.0 + 0.3345167463*(pow(x,2.0) +  pow(y,2.0)) +
              0.09871707060*(pow(x,4.0) +  pow(y,4.0))),2.0));
}





