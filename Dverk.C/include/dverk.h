/*
 *  dverk.h
 *  Dverk-Test
 *
 *  Created by ashley on 15/05/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

/*
c***********************************************************************
c                                                                      *
c  summary of the components of the communications vector              *
c                                                                      *
c     prescribed at the option       determined by the program         *
c           of the user                                                *
c                                                                      *
c                                    c(10) rreb(rel roundoff err bnd)  *
c     c(1) error control indicator   c(11) dwarf (very small mach no)  *
c     c(2) floor value               c(12) weighted norm y             *
c     c(3) hmin specification        c(13) hmin                        *
c     c(4) hstart specification      c(14) hmag                        *
c     c(5) scale specification       c(15) scale                       *
c     c(6) hmax specification        c(16) hmax                        *
c     c(7) max no of fcn evals       c(17) xtrial                      *
c     c(8) interrupt no 1            c(18) htrial                      *
c     c(9) interrupt no 2            c(19) est                         *
c                                    c(20) previous xend               *
c                                    c(21) flag for xend               *
c                                    c(22) no of successful steps      *
c                                    c(23) no of successive failures   *
c                                    c(24) no of fcn evals             *
c                                                                      *
c  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        */
#include"control.h"
struct _communication {
        int n;
        
        control *ctrl;
        void (*solver)();
        
//       prescribed at the option of the user
        int errorcontrol;                                                                                                       //      c[0]    error control indicator
        double floorvalue;                                                                                                      //      c[1]    floor value
        double hminspec;                                                                                                        //      c[2]    hmin specification
        double hstart;                                                                                                          //      c[3]    hstart specification
        double scalespec;                                                                                                       //      c[4]    scale specification
        double hmaxspec;                                                                                                        //      c[5]    hmax specification
        int maxevals;                                                                                                           //      c[6]    max no of fcn evals
        int interrupt1;                                                                                                         //      c[7]    interrupt no 1
        int interrupt2;                                                                                                         //      c[8]    interrupt no 2

//      determined by the programme
        double rreb;                                                                                                            //      c[9]    rreb(rel roundoff err bnd)
        double dwarf;                                                                                                           //      c[10]   dwarf (very small mach no)
        double weightednormy;                                                                                           //      c[11]   weighted norm y
        double hmin;                                                                                                            //      c[12]   hmin
        double hmag;                                                                                                            //      c[13]   hmag
        double scale;                                                                                                           //      c[14]   scale
        double hmax;                                                                                                            //      c[15]   hmax
        double xtrial;                                                                                                          //      c[16]   xtrial
        double htrial;                                                                                                          //      c[17]   htrial
        double est;                                                                                                                     //      c[18]   est
        double prevxend;                                                                                                        //      c[19]   previous xend
        int xendflag;                                                                                                           //      c[20]   flag for xend
        int nsuccessful;                                                                                                        //      c[21]   no of successful steps
        int nfailed;                                                                                                            //      c[22]   no of successive failures
        int nevaluations;                                                                                                       //      c[23]   no of fcn evals
        double *floorvalues;                                                                                            //      c[30..n+29]     if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values
};

typedef struct _communication dverkcomm;

dverkcomm *makeDverkCommunication(control *c);
void freeDverkCommunication(dverkcomm *comm);

//void discrete(dverkcomm *c, void (*fcn)(), double x, double *y, double xend, double tol, int *ind, int nw, void *fdata);
//void dverk(dverkcomm *c, void (*fcn)(), double x, double *y, double xend, double tol, int *ind, int nw, void *fdata);