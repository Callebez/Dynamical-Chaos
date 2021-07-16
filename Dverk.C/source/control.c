/*
 *  control.c
 *  dynamo.cgi
 *
 *  Created by ashley on 19/05/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/control.h"

// void output(int n, double *y, double t);

char *astrcpy(char *source)
{
        int n = strlen(source);
        char *dest = malloc(n + 1);
        strcpy(dest, source);
        
        return dest;
}

constant *makeconstant()
{
        constant *con = (constant *)malloc(sizeof(constant));
        con->name = NULL;
        con->value = 0.0;
        con->info = NULL;
        con->row = -1;
        con->column = -1;
        
        return con;
}

void freeconstant(constant *con)
{
        if (con->name!=NULL) free(con->name);
        if (con->info!=NULL) free(con->info);
        
        free(con);
}

variable *makevariable()
{
        variable *var = (variable *)malloc(sizeof(variable));
        var->name = NULL;
        var->linecolour = NULL;
        var->markercolour = NULL;
        
        return var;
}

void freevariable(variable *var)
{
        free(var->name);
        if (var->linecolour!=NULL) free(var->linecolour);
        if (var->markercolour!=NULL) free(var->markercolour);
        free(var);
}

plot* makeplot(plottype t, int nlines, int ordinate, int *lines)
{
        plot *p = (plot *)calloc(1, sizeof(plot));
        
        p->type = t;
        p->nlines = nlines;
        p->ordinate = ordinate;
        p->abscissae = (int *)calloc(nlines, sizeof(int));
        int i;
        for(i=0; i<nlines; i++) {
                p->abscissae[i] = lines[i];
        }
        
        p->label = "";
        p->title = "";
        p->xlabel = "";
        p->ylabel = "";
        
        return p;
}

void freeplot(plot *p)
{
        free(p->abscissae);
        free(p);
}

control *makecontrol(int n, int p, int m)
{
        control *c = (control *)calloc(1, sizeof(control));

        c->title = "";
        c->description = "";

        c->output = output;
        c->outputcounter = 0;
        c->maxoutput = 0;
        c->results = NULL;

        c->nparameters = p;
        c->nvariables = n;
        
//      c->variablenames = (char **)malloc(n*sizeof(char *));
//      c->variablevalues = (double *)calloc(n, sizeof(double));
        c->variables = (variable **)malloc(n*sizeof(variable *));
        int i;
        for(i=0; i<n; i++) {
                c->variables[i] = makevariable();
        }
        
/*      c->parameternames = (char **)malloc(p*sizeof(char *));
        c->parameterinfo = (char **)malloc(p*sizeof(char *));
        c->parametervalues = (double *)calloc(p, sizeof(double));
        for(i=0; i<p; i++) {
                c->parameternames[i] = "";
                c->parameterinfo[i] = "";
        } */
        c->parameters = (constant **)malloc(p*sizeof(constant *));
        for(i=0; i<p; i++) {
                c->parameters[i] = makeconstant();
                c->parameters[i]->column = 0;
                c->parameters[i]->row = i;
        }
        
        c->nplots = m;
        c->plots = (plot **)malloc(p*sizeof(plot *));
        for(i=0; i<m; i++) {
                c->plots[i] = NULL;
        }
        
        c->ncolumns = 1;
        c->nrows = c->nparameters;
        
        return c;
}

void freecontrol(control *c)
{
        if (c->results!=NULL) {
                free(c->results);
        }

//      free(c->variablevalues);
//      free(c->parametervalues);
        int i;
        for(i=0; i<c->nvariables; i++) {
//              free(c->variablenames[i]);
                freevariable(c->variables[i]);
        }
        free(c->variables);
        
        for(i=0; i<c->nparameters; i++) {
//              free(c->parameternames[i]);
//              free(c->parameterinfo[i]);
                freeconstant(c->parameters[i]);
        }
        free(c->parameters);
        
        for(i=0; i<c->nplots; i++) {
                if(c->plots[i] != NULL) freeplot(c->plots[i]);
        }
        free(c->plots);
        free(c);
}

void initvariables(control *c, double *y)
{
        int i;
        for(i=0; i<c->nvariables; i++) {
                y[i] = c->variables[i]->value;
        }
}

void output(control *c, int n, double *y, double t)
{
//      printf("%g", t);
        
        int p = c->outputcounter*(c->nvariables + 1);
        c->results[p++] = t;
        int i;
        for(i=0; i<n; i++) {
//              printf("\t%g", y[i]);
                c->results[p++] = y[i];
        }
//      printf("\n");
        
        
        (c->outputcounter)++;

}

void printresults(control *c)
{
        int i;
        int j;
        
        for(i=0; i<c->outputcounter; i++) {
                for(j=0; j<=c->nvariables; j++) {
                        printf("%g ", c->results[i*(c->nvariables + 1) + j]);
                }
                printf("\n");
        }
}
