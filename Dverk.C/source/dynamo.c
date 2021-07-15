/*
 *  dynamo.c
 *  dynamo.cgi
 *
 *  Created by ashley on 18/05/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "../include/control.h"
#include "../include/dynamo.h"

//void initialise(int *nconstants, int *nvariables, int *nplots);
// void initcontrol(control *c);
void fcndynamo(int n, double t, double *y, double *yprime, void *fdata);

/*void initialise(int *nconstants, int *nvariables, int *nplots)
{
        *nconstants = 4;
        *nvariables = 3;
        *nplots = 4;
} */

control *makedynamo()
{
        int nconstants = 4;
        int nvariables = 3;
        int nplots = 4;
        control *c = makecontrol(nvariables, nconstants, nplots);
        
        c->title = "Double Disc Dynamo ODE model";
        c->description = "Adapted from the gnu graph 'ode' package.";
        
        c->integrator = dverkintegrator;
//      c->init = initcontrol;
        c->fcn = fcndynamo;
        
        c->parameters[0]->value = 1.0;
        c->parameters[0]->name = astrcpy("A");
        c->parameters[0]->info = astrcpy("A parameter");
        
        c->parameters[1]->value = 1.0;
        c->parameters[1]->name = astrcpy("V");
        c->parameters[1]->info = astrcpy("V parameter");
        
        c->parameters[2]->value = 14.625;
        c->parameters[2]->name = astrcpy("Q");
        c->parameters[2]->info = astrcpy("Q parameter");
        
        c->parameters[3]->value = 5.0;
        c->parameters[3]->name = astrcpy("S");
        c->parameters[3]->info = astrcpy("S parameter");

        c->variables[0]->name = astrcpy("w");
        c->variables[1]->name = astrcpy("y");
        c->variables[2]->name = astrcpy("z");
        c->variables[0]->value = 1.0;
        c->variables[1]->value = 1.0;
        c->variables[2]->value = 1.0;
        
        int lines[] = {1, 2, 3};
        c->plots[0] = makeplot(lineplot, 3, 0, lines);
        c->plots[0]->title = astrcpy("Time series");
        c->plots[0]->label = astrcpy("w-y-z Time Series Plot");
        c->plots[0]->xlabel = astrcpy("y");
        c->plots[0]->ylabel = astrcpy("t");

        int lines2[] = {3};
        c->plots[1] = makeplot(lineplot, 1, 2, lines2);
        c->plots[1]->title = astrcpy("Phase Space");
        c->plots[1]->label = astrcpy("z-y Phase Space Plot");
        c->plots[1]->xlabel = astrcpy("w");
        c->plots[1]->ylabel = astrcpy("y");

        int lines3[] = {2};
        c->plots[2] = makeplot(lineplot, 1, 1, lines3);
        c->plots[2]->title = astrcpy("Phase Space");
        c->plots[2]->label = astrcpy("w-y Phase Space Plot");
        c->plots[2]->xlabel = astrcpy("w");
        c->plots[2]->ylabel = astrcpy("y");
        
        int lines4[] = {3};
        c->plots[3] = makeplot(lineplot, 1, 1, lines4);
        c->plots[3]->title = astrcpy("Phase Space");
        c->plots[3]->label = astrcpy("z-w Phase Space Plot");
        c->plots[3]->xlabel = astrcpy("w");
        c->plots[3]->ylabel = astrcpy("y");
        
        c->t0 = 0.0;
        c->t1 = 25.0;
        c->dt = 0.05;
        c->tol = 1.0e-07;
        
        return c;
}

/*void initcontrol(control *c)
{
        c->title = "Double Disc Dynamo ODE model";
        c->description = "Adapted from the gnu graph 'ode' package.";
        
        c->integrator = dverkintegrator;

        c->parametervalues[0] = 1.0;
        c->parameternames[0] = astrcpy("A");
        c->parameterinfo[0] = astrcpy("A parameter");
        
        c->parametervalues[1] = 1.0;
        c->parameternames[1] = astrcpy("V");
        c->parameterinfo[1] = astrcpy("V parameter");
        
        c->parametervalues[2] = 14.625;
        c->parameternames[2] = astrcpy("Q");
        c->parameterinfo[2] = astrcpy("Q parameter");
        
        c->parametervalues[3] = 5.0;
        c->parameternames[3] = astrcpy("S");
        c->parameterinfo[3] = astrcpy("S parameter");

        c->variablenames[0] = astrcpy("w");
        c->variablenames[1] = astrcpy("y");
        c->variablenames[2] = astrcpy("z");
        c->variablevalues[0] = 1.0;
        c->variablevalues[1] = 1.0;
        c->variablevalues[2] = 1.0;
        
        int lines[] = {1, 2, 3};
        c->plots[0] = makeplot(lineplot, 3, 0, lines);
        c->plots[0]->title = astrcpy("Time series");
        c->plots[0]->label = astrcpy("w-y-z Time Series Plot");
        c->plots[0]->xlabel = astrcpy("y");
        c->plots[0]->ylabel = astrcpy("t");

        int lines2[] = {3};
        c->plots[1] = makeplot(lineplot, 1, 2, lines2);
        c->plots[1]->title = astrcpy("Phase Space");
        c->plots[1]->label = astrcpy("z-y Phase Space Plot");
        c->plots[1]->xlabel = astrcpy("w");
        c->plots[1]->ylabel = astrcpy("y");

        int lines3[] = {2};
        c->plots[2] = makeplot(lineplot, 1, 1, lines3);
        c->plots[2]->title = astrcpy("Phase Space");
        c->plots[2]->label = astrcpy("w-y Phase Space Plot");
        c->plots[2]->xlabel = astrcpy("w");
        c->plots[2]->ylabel = astrcpy("y");
        
        int lines4[] = {3};
        c->plots[3] = makeplot(lineplot, 1, 1, lines4);
        c->plots[3]->title = astrcpy("Phase Space");
        c->plots[3]->label = astrcpy("z-w Phase Space Plot");
        c->plots[3]->xlabel = astrcpy("w");
        c->plots[3]->ylabel = astrcpy("y");
        
        c->t0 = 0.0;
        c->t1 = 25.0;
        c->dt = 0.05;
        c->tol = 1.0e-07;
} */

void fcndynamo(int n, double t, double *y, double *yprime, void *fdata)
{
//      yprime[0] = -y[1] - 0.0*y[0];
//      yprime[1] = y[0];
//      yprime[0] = -y[0];
        control *c = (control *)fdata;
        constant **params = c->parameters;
        double A = params[0]->value;
        double V = params[1]->value;
        double Q = params[2]->value;
        double S = params[3]->value;

        yprime[0] = Q - y[2] * y[1] - V * y[0];                                                                 // w = Q - z * y - V * w
        yprime[1] = S * ( A * y[2] - y[1]);                                                                             // y = S * ( A * z - y)
        yprime[2] = y[0] * y[1] - y[2];                                                                                 // z = w * y - z
}
