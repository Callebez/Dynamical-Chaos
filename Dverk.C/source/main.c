#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "../include/cgic.h"
#include "../include/control.h"
#include "../include/dverk.h"

#include "../include/dynamo.h"
#include "../include/rossler.h"
#include "../include/2sp-osc.h"
#include "../include/rinaldi.h"

#define STRMAX 256

/*void initialise(int *nconstants, int *nvariables, int *nplots);
void initcontrol(control *c);
void fcn(int n, double t, double *y, double *yprime, void *fdata);
void output(int n, double *y, double t);
control *makeproblem(); */

int cgiMain () 
{
        int i;

        cgiHeaderContentType("text/json");
//      printf("Header line (ignored)\n");
        printf("{\"header\":\"odemodels.cgi\"\n");

        char *errorstring = NULL;

        char modelName[STRMAX];
        control *c = NULL;
        int nProbems = 4;
        char *problemList[] = {"dynamo", "rossler", "2sp-osc", "rinaldi"};
        control *(*makeFunctions[])() = {makedynamo, makerossler, make2sposc, makerinaldi};
        
        cgiFormStringNoNewlines("problem", modelName, STRMAX);
        
        int problemNumber = -1;
        for(i=0; i<nProbems; i++) {
                if(strcmp(problemList[i], modelName)==0) {
                        problemNumber = i;
                        break;
                }
        }
        
/*      if (strncmp(modelName, "dynamo", STRMAX)==0) {
                c = makedynamo();
        } else if (strncmp(modelName, "rossler", STRMAX)==0) {
                c = makerossler();
        } else if (strncmp(modelName, "2sp-osc", STRMAX)==0) {
                c = make2sposc();
        } else if (strncmp(modelName, "rinaldi", STRMAX)==0) {
                c = makerinaldi(); */
        if (problemNumber>=0) {
                c = makeFunctions[problemNumber]();
        } else {
                printf(",\"problems\":[");
                for(i=0; i<nProbems; i++) {
                        if (i>0) {
                                printf(", ");
                        }
                        printf("\"%s\"", problemList[i]);
                }
                printf("]\n");
//              errorstring = astrcpy("No/invalid problem specified.");
        }
        
        if (c!=NULL) {
                cgiFormDoubleBounded("starttime", &(c->t0), 0.0, 1.0e38, c->t0);
                cgiFormDoubleBounded("stoptime", &(c->t1), 0.0, 1.0e38, c->t1);
                cgiFormDoubleBounded("timestep", &(c->dt), 0.0, 1.0e38, c->dt);
                cgiFormDoubleBounded("tolerance", &(c->tol), 0.0, 1.0e38, c->tol);
                
                for(i=0; i<c->nparameters; i++) {
                        cgiFormDoubleBounded(c->parameters[i]->name, &(c->parameters[i]->value), -1.0e38, 1.0e38, c->parameters[i]->value);
                        //              fprintf(cgiOut, "c[%d] = %f<BR>\n", i, out->c[i]);
                }
        }
        
        cgiFormStringNoNewlines("header", modelName, STRMAX);
        int printheader = strncmp(modelName, "yes", STRMAX)==0;
        if (printheader && (c!=NULL)) {
//      Header
                printf(",\"modelname\":\"%s\",\n", c->title);
//              printf("%d\n", 1);
                printf("\"description\":[\"%s\"],\n", c->description);
                printf("\"integratortype\":%d,\n", (int)c->integrator);
                printf("\"controlparams\":[%g, %g, %g, %g],\n", c->t0, c->t1, c->dt, c->tol);
                printf("\"nrows\":%d,\n", c->nrows);
                printf("\"ncolumns\":%d,\n", c->ncolumns);
                
//      Constants
                printf("\"nvars\":%d,\n", c->nvariables);
                printf("\"nconsts\":%d,\n", c->nparameters);
                printf("\"constants\":[\n");
                for(i=0; i<c->nparameters; i++) {
                        if (i>0) {
                                printf(",");
                        }
                        printf("{\"val\":%g, \n\"name\":\"%s\", \n\"info\":\"%s\",\n", c->parameters[i]->value, c->parameters[i]->name, c->parameters[i]->info);
                        printf("\"row\":%d, ", c->parameters[i]->row);
                        printf("\"column\":%d", c->parameters[i]->column);
                        printf("}\n");
                }
                printf("],\n");
                
//      Variables
/*              printf("%g\n%s\n%s\n", y[0], "w", "w");
                printf("%g\n%s\n%s\n", y[0], "y", "y");
                printf("%g\n%s\n%s\n", y[0], "z", "z"); */
                char *colours[] = {"black", "blue", "red", "green", "yellow", "pink", "salmon", "aqua"};
                int ncolour = 0;
                int ncolours = 8;
                printf("\"variables\":[{\"name\":\"t\",\n \"lineColour\": \"black\",\n \"lineStyle\": 0,\n \"markerColour\":\"\",\n");
                printf(" \"markerStyle\":0,\n \"markerSize\":0,\n \"lineThickness\":1}\n");
                for(i=0; i<c->nvariables; i++) {
                        printf(",{");
                        char *colour = colours[ncolour++];
                        printf("\"name\":\"%s\",\n", c->variables[i]->name);
                        printf(" \"lineColour\": \"%s\",\n", colour);
                        printf(" \"lineStyle\": %d,\n", 1);
                        printf(" \"markerColour\":\"%s\",\n", colour);
                        printf(" \"markerStyle\":0,\n");
                        printf(" \"markerSize\":%d,\n", 5);
                        printf(" \"lineThickness\":%d", 1);
                        printf("}\n");
                        if (ncolour==ncolours) ncolour = 0;
                }
                printf("],\n");

//      Plots
                printf("\"nplots\":%d,\n\"plots\":[\n", c->nplots);
                for(i=0; i<c->nplots; i++) {
                        if (i>0) {
                                printf(",");
                        }
                        printf("{\"title\":\"%s\",\n", c->plots[i]->title);                                     // Title
                        printf("\"label\":\"%s\",\n", c->plots[i]->label);                                      // Label for selector
                        printf("\"xlabel\":\"%s\",\n", c->plots[i]->ylabel);                            //      y-label
                        printf("\"ylabel\":\"%s\",\n", c->plots[i]->xlabel);                            // x-label
//                      printf("\"%d\",\n", c->plots[i]->nlines);
                        printf("\"lines\":[%d,", c->plots[i]->ordinate);
                        int j;
                        for(j=0; j<c->plots[i]->nlines; j++) {
                                if (j>0) {
                                        printf(",");
                                }
                                printf(" %d", c->plots[i]->abscissae[j]);               // variables to print (0 is time)
                        }
                        printf("]}\n");
                }
                printf("]\n");
                
//      } else {
        }

        dverkcomm *dvc = NULL;
        double *y = NULL;
        if (c!=NULL) {
                y = calloc(c->nvariables, sizeof(double));
                initvariables(c, y);
                dvc = makeDverkCommunication(c);
                dvc->ctrl = c;
                
                c->maxoutput = floor(1.5*(c->t1 - c->t0)/c->dt);
                c->results = (double *)calloc((c->nvariables + 1)*c->maxoutput, sizeof(double));
                c->outputcounter = 0;
                
//              double t = c->t0;
//              double tol = c->tol;
//              int n = c->nvariables;
                c->output(c, c->nvariables, y, c->t0);
        }

/*      while(t<=c->t1) {
                double tend = t + c->dt;
                dvc->solver(dvc, c->fcn, t, y, tend, tol, &ind, n, c);
                t = tend;
                if (t<=c->t1) c->output(n, y, t);
        } */
        
//      printresults(c);
//      int comma = 0;
        
        cgiFormStringNoNewlines("results", modelName, STRMAX);
        if ((strncmp(modelName, "yes", STRMAX)==0) && (c!=NULL)) {
                int ind = 1;
                dvc->solver(dvc, c->fcn, c->t0, y, c->t1, c->tol, &ind, c->nvariables, c);
                printf(",\"results\":[");
                int j;
                for(j=0; j<=c->nvariables; j++) {
                        double minimum;
                        double maximum;
                        printf("{");
                        printf("\"data\":[");
                        for(i=0; i<c->outputcounter; i++) {
                                double y = c->results[i*(c->nvariables + 1) + j];
                                if (i>0) {
                                        printf(", ");
                                        if (y>maximum) {
                                                maximum = y;
                                        }
                                        if (y<minimum) {
                                                minimum = y;
                                        }
                                } else {
                                        maximum = y;
                                        minimum = y;
                                }
                                printf("%g", y);
                        }
                        printf("]");
                        printf(", \"minimum\": %g", minimum);
                        printf(", \"maximum\": %g", maximum);
                        printf(", \"increment\": %g", (maximum - minimum)/4.0);
                        printf("}");
                        if (j<c->nvariables) {
                                printf(", ");
                        }
                        printf("\n"); 
                }
                printf("]\n");
        }
        printf(",\"errors\":[");
        if (errorstring!=NULL) printf("\"%s\"", errorstring);
        printf("]\n");
        printf("}\n");

        if (dvc!=NULL) freeDverkCommunication(dvc);
        if (y!=NULL) free(y);
        if (errorstring!=NULL) free(errorstring);

    return 0;
}