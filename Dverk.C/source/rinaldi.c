/*
 *  rinaldi.c
 *  odemodels.cgi
 *
 *  Created by ashley on 09/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "../include/control.h"
#include "../include/rinaldi.h"

void fcnrinaldi(int n, double t, double *y, double *yprime, void *fdata);

control *makerinaldi()
{
        int nconstants = 13;
        int nvariables = 3;
        int nplots = 4;
        control *c = makecontrol(nvariables, nconstants, nplots);
        
        c->title = "Love Dynamics";
        c->description = "Rinaldi's model of the love between Petrarch & Laura";

        c->integrator = dverkintegrator;
        c->fcn = fcnrinaldi;

        int p = 0;
        c->parameters[p]->value = 3.0;
        c->parameters[p]->name = astrcpy("alpha1");
        c->parameters[p]->info = astrcpy("alpha1 parameter");
        p++;
        
        c->parameters[p]->value = 1.0;
        c->parameters[p]->name = astrcpy("alpha2");
        c->parameters[p]->info = astrcpy("alpha2 parameter");
        p++;
        
        c->parameters[p]->value = 0.1;
        c->parameters[p]->name = astrcpy("alpha3");
        c->parameters[p]->info = astrcpy("alpha3 parameter");
        p++;
        
        c->parameters[p]->value = 1.0;
        c->parameters[p]->name = astrcpy("beta1");
        c->parameters[p]->info = astrcpy("beta1 parameter");
        p++;
        
        c->parameters[p]->value = 5.0;
        c->parameters[p]->name = astrcpy("beta2");
        c->parameters[p]->info = astrcpy("beta2 parameter");
        p++;
        
        c->parameters[p]->value = 0.1e+02;
        c->parameters[p]->name = astrcpy("beta3");
        c->parameters[p]->info = astrcpy("beta3 parameter");
        p++;
        
        c->parameters[p]->value = 1.0;
        c->parameters[p]->name = astrcpy("gamma");
        c->parameters[p]->info = astrcpy("gamma parameter");
        p++;
        
        c->parameters[p]->value = 1.0;
        c->parameters[p]->name = astrcpy("delta");
        c->parameters[p]->info = astrcpy("delta parameter");
        p++;
        
        c->parameters[p]->value = 2.0;
        c->parameters[p]->name = astrcpy("A_L");
        c->parameters[p]->info = astrcpy("A_L parameter");
        p++;
        
        c->parameters[p]->value = -1.0;
        c->parameters[p]->name = astrcpy("A_P");
        c->parameters[p]->info = astrcpy("A_P parameter");
        p++;
        
        c->parameters[p]->value = 0.0;
        c->parameters[p]->name = astrcpy("L(0)");
        c->parameters[p]->info = astrcpy("Initial L");
        p++;
        
        c->parameters[p]->value = 0.0;
        c->parameters[p]->name = astrcpy("P(0)");
        c->parameters[p]->info = astrcpy("Initial P");
        p++;
        
        c->parameters[p]->value = 0.0;
        c->parameters[p]->name = astrcpy("Z(0)");
        c->parameters[p]->info = astrcpy("Initial Z");
        p++;
        
        c->variables[0]->name = astrcpy("L");
        c->variables[1]->name = astrcpy("P");
        c->variables[2]->name = astrcpy("Z");
        c->variables[0]->value = c->parameters[10]->value;
        c->variables[1]->value = c->parameters[11]->value;
        c->variables[2]->value = c->parameters[12]->value;
        
        int lines[] = {1, 2, 3};
        c->plots[0] = makeplot(lineplot, 3, 0, lines);
        c->plots[0]->title = astrcpy("Time series");
        c->plots[0]->label = astrcpy("x-y-z Time Series Plot");
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
        c->plots[2]->label = astrcpy("x-y Phase Space Plot");
        c->plots[2]->xlabel = astrcpy("x");
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

void fcnrinaldi(int n, double t, double *y, double *yprime, void *fdata)
{
        control *c = (control *)fdata;
        constant **params = c->parameters;
        double alpha1 = params[0]->value;
        double alpha2 = params[1]->value;
        double alpha3 = params[2]->value;

        double beta1 = params[3]->value;
        double beta2 = params[4]->value;
        double beta3 = params[5]->value;

        double gamma = params[6]->value;
        double delta = params[7]->value;

        double A_L = params[8]->value;
        double A_P = params[9]->value;

        double L = y[0];
        double P = y[1];
        double Z = y[2];
        
        double pp = P/gamma;
        yprime[0] = -alpha1*L + beta1*(P*(1 - pp*pp) + A_P);
        yprime[1] = -alpha2*P + beta2*(L + A_L/(1 + delta*Z));
        yprime[2] = -alpha3*Z + beta3*P;
}
