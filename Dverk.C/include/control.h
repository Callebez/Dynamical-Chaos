/*
 *  control.h
 *  dynamo.cgi
 *
 *  Created by ashley on 19/05/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

enum _plottype {lineplot};
typedef enum _plottype plottype;

enum _integratortype {discreteintegrator, eulerintegrator, dde23integrator, dverkintegrator};
typedef enum _integratortype integratortype;

struct _plot {
        plottype type;
        int nlines;
        int ordinate;
        int *abscissae;
        char *label;
        char *title;
        char *xlabel;
        char *ylabel;
};
typedef struct _plot plot;

struct _constant {
        char *name;
        double value;
        char *info;
        int row;
        int column;
};
typedef struct _constant constant;

struct _variable {
        char *name;
        double value;
        
        char *linecolour;
        int lineStyle;
        char *markercolour;
        int markerstyle;
        int markersize;
        int linethickness;
};
typedef struct _variable variable;

struct _control {
        char *title;
        char *description;
        
        integratortype integrator;
        void (*fcn)();
        void (*init)();
        void (*output)();

        int nparameters;
        int nvariables;
//      char **variablenames;
//      double *variablevalues;
        variable **variables;
/*      char **parameternames;
        char **parameterinfo;
        double *parametervalues; */
        constant **parameters;
        
        int outputcounter;
        int maxoutput;
        double *results;
        
        int ncolumns;
        int nrows;
        
        int nplots;
        plot **plots;
        
        double t0;
        double t1;
        double dt;
        double tol;
};
typedef struct _control control;

char *astrcpy(char *source);
constant *makeconstant();
void freeconstant(constant *con);
variable *makevariable();
void freevariable(variable *var);
plot* makeplot(plottype t, int nlines, int ordinate, int *lines);
control *makecontrol(int n, int p, int m);
void freecontrol(control *c);
void initvariables(control *c, double *y);
void output(control *c, int n, double *y, double t);
void printresults(control *c);