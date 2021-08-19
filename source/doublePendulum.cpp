#include "../include/doublePendulum.hpp"
std::vector<double> doublePendulum(std::vector<double> coord, double inertia)
{
    std::vector<double> xcoord (4,0);

    xcoord[0] = (6.0/inertia)*(2.0*coord[2]-3.0*cos(coord[0]-coord[1])*coord[3])/(16.0-9.0*cos(coord[0]-coord[1]));
    xcoord[1] = (6.0/inertia)*(8.0*coord[2]-3.0*cos(coord[0]-coord[1])*coord[3])/(16.0-9.0*cos(coord[0]-coord[1]));
    xcoord[2] = -0.5*inertia*(xcoord[0]*xcoord[1]*sin(coord[0]-coord[1]) + 3 * (9.81/inertia)*sin(coord[0]) );
    xcoord[3] = -0.5*inertia*(-xcoord[0]*xcoord[1]*sin(coord[0]-coord[1]) +  (9.81/inertia)*sin(coord[1]) );


    return xcoord;
}