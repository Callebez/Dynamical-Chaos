/*
 *  2sp-osc.c
 *  dynamo.cgi
 *
 *  Created by ashley on 22/08/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

//#include "2sp-osc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "../include/control.h"
#include "../include/2sp-osc.h"

void fcn2sposc(int n, double t, double *y, double *yprime, void *fdata);

double _auxiliary0;

control *make2sposc()
{
        int nconstants = 2;
        int nvariables = 2;
        int nplots = 2;
        control *c = makecontrol(nvariables, nconstants, nplots);
        
        c->title = "Two Species Oscillator";
        c->description = "Two Species Biochemical Oscillator";

        c->integrator = dverkintegrator;
        c->fcn = fcn2sposc;

    c->variables[0]->value = 1;
    c->variables[1]->value = 1;
    c->parameters[0]->value = 0.1;
    c->parameters[1]->value = 0.5;

    c->parameters[0]->name = astrcpy("a");
    c->parameters[0]->info = astrcpy("Parameter 0");
    c->parameters[1]->name = astrcpy("b");
    c->parameters[1]->info = astrcpy("Parameter 1");

    c->variables[0]->name = astrcpy("x");
    c->variables[1]->name = astrcpy("y");

        int lines[] = {1, 2};
        c->plots[0] = makeplot(lineplot, 2, 0, lines);
        c->plots[0]->title = astrcpy("Time series");
        c->plots[0]->label = astrcpy("x-y Time Series Plot");
        c->plots[0]->xlabel = astrcpy("x");
        c->plots[0]->ylabel = astrcpy("t");

        int lines2[] = {2};
        c->plots[1] = makeplot(lineplot, 1, 1, lines2);
        c->plots[1]->title = astrcpy("Phase Space");
        c->plots[1]->label = astrcpy("x-y Phase Space Plot");
        c->plots[1]->xlabel = astrcpy("x");
        c->plots[1]->ylabel = astrcpy("y");

    c->t0 = 0.0;
    c->t1 = 100;
    c->dt = 1.0;
    c->tol = 1.0e-08;
        
        return c;
}

void fcn2sposc(int n, double t, double *s, double *g, void *fdata)
{
    control *ctrl = (control *)fdata;
    constant **c = ctrl->parameters;
    _auxiliary0 = (s[0]*s[0])*s[1];
    g[0] = (c[0]->value - s[0])+_auxiliary0;
    g[1] = c[1]->value - _auxiliary0;
}