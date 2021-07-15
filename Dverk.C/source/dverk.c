/*
 *  dverk.c
 *  Dverk-Test
 *
 *  Created by ashley on 15/05/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "../inculde/control.h"
#include "../inculde/dverk.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void discrete(dverkcomm *c, void (*fcn)(), double x, double *y, double xend, double tol, int *ind, int nw, void *fdata);
void dverk(dverkcomm *c, void (*fcn)(), double x, double *y, double xend, double tol, int *ind, int nw, void *fdata);

double dmax1(double x, double y)
{
        double z = (x>y)?x:y;
        return z;
}

double dmin1(double x, double y)
{
        double z = (x<y)?x:y;
        return z;
}

double dsign(double x, double y)
{
        double z;
        if (y>0) {
                z = fabs(x);
        } else {
                z = -fabs(x);
        }
        return z;
}

dverkcomm *makeDverkCommunication(control *c)
{
        int n = c->nvariables;
        dverkcomm *comm = (dverkcomm *)calloc(1, sizeof(dverkcomm));
        comm->n = n;
        if (c->integrator==discreteintegrator) {
                comm->solver = discrete;
        } else if (c->integrator==dverkintegrator) {
                comm->floorvalues = calloc(n, sizeof(double));
                comm->solver = dverk;
        } else {
        }
        return comm;
}

void freeDverkCommunication(dverkcomm *comm)
{
        free(comm->floorvalues);
        freecontrol(comm->ctrl);
        free(comm);
}

void discrete(dverkcomm *c, void (*fcn)(), double x, double *y, double xend, double tol, int *ind, int nw, void *fdata)
{
}

void euler(dverkcomm *c, void (*fcn)(), double x, double *y, double xend, double tol, int *ind, int nw, void *fdata)
{
        double t = x;
        int n = c->n;
        double tout = x + c->ctrl->dt;
        double tstep = c->hstart;
        double *w0 = (double *)calloc(n, sizeof(double));

        c->ctrl->output(c->ctrl, n, y, t);

        while (t<=xend) {
                fcn(n, x, y, w0, fdata);
                int j;
                for(j=0; j<n; j++) {
                        y[j] += tstep*w0[j];
                }
                t += tstep;

                while ((x>tout) && (c->ctrl->dt>0.0)) {
                        c->ctrl->output(c->ctrl, n, y, x);
                        tout += c->ctrl->dt;
                }
        }
        
        free(w0);
}

//      subroutine dverk (n, fcn, x, y, xend, tol, ind, c, nw, w)
//      integer n, ind, nw, k
//      double precision x, y(n), xend, tol, c(24), w(nw,9), temp

/*
c***********************************************************************
c                                                                      *
c note added 11/14/85.                                                 *
c                                                                      *
c if you discover any errors in this subroutine, please contact        *
c                                                                      *
c        kenneth r. jackson                                            *
c        department of computer science                                *
c        university of toronto                                         *
c        toronto, ontario,                                             *
c        canada   m5s 1a4                                              *
c                                                                      *
c        phone: 416-978-7075                                           *
c                                                                      *
c        electronic mail:                                              *
c        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
c        csnet:  krj@toronto                                           *
c        arpa:   krj.toronto@csnet-relay                               *
c        bitnet: krj%toronto@csnet-relay.arpa                          *
c                                                                      *
c dverk is written in fortran 66.                                      *
c                                                                      *
c the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
c set for a  vax  in  double  precision.  they  should  be  reset,  as *
c described below, if this program is run on another machine.          *
c                                                                      *
c the c array is declared in this subroutine to have one element only, *
c although  more  elements  are  referenced  in this subroutine.  this *
c causes some compilers to issue warning messages.  there is,  though, *
c no  error  provided  c is declared sufficiently large in the calling *
c program, as described below.                                         *
c                                                                      *
c the following external statement  for  fcn  was  added  to  avoid  a *
c warning  message  from  the  unix  f77 compiler.  the original dverk *
c comments and code follow it.                                         *
c                                                                      *
c***********************************************************************
c
      external fcn
c
c***********************************************************************
c                                                                      *
c     purpose - this is a runge-kutta  subroutine  based  on  verner's *
c fifth and sixth order pair of formulas for finding approximations to *
c the solution of  a  system  of  first  order  ordinary  differential *
c equations  with  initial  conditions. it attempts to keep the global *
c error proportional to  a  tolerance  specified  by  the  user.  (the *
c proportionality  depends  on the kind of error control that is used, *
c as well as the differential equation and the range of integration.)  *
c                                                                      *
c     various options are available to the user,  including  different *
c kinds  of  error control, restrictions on step sizes, and interrupts *
c which permit the user to examine the state of the  calculation  (and *
c perhaps make modifications) during intermediate stages.              *
c                                                                      *
c     the program is efficient for non-stiff systems.  however, a good *
c variable-order-adams  method  will probably be more efficient if the *
c function evaluations are very costly.  such a method would  also  be *
c more suitable if one wanted to obtain a large number of intermediate *
c solution values by interpolation, as might be the case  for  example *
c with graphical output.                                               *
c                                                                      *
c                                    hull-enright-jackson   1/10/76    *
c                                                                      *
c***********************************************************************
c                                                                      *
c     use - the user must specify each of the following                *
c                                                                      *
c     n  number of equations                                           *
c                                                                      *
c   fcn  name of subroutine for evaluating functions - the  subroutine *
c           itself must also be provided by the user - it should be of *
c           the following form                                         *
c              subroutine fcn(n, x, y, yprime)                         *
c              integer n                                               *
c              double precision x, y(n), yprime(n)                     *
c                      *** etc ***                                     *
c           and it should evaluate yprime, given n, x and y            *
c                                                                      *
c     x  independent variable - initial value supplied by user         *
c                                                                      *
c     y  dependent variable - initial values of components y(1), y(2), *
c           ..., y(n) supplied by user                                 *
c                                                                      *
c  xend  value of x to which integration is to be carried out - it may *
c           be less than the initial value of x                        *
c                                                                      *
c   tol  tolerance - the subroutine attempts to control a norm of  the *
c           local  error  in  such  a  way  that  the  global error is *
c           proportional to tol. in some problems there will be enough *
c           damping  of  errors, as well as some cancellation, so that *
c           the global error will be less than tol. alternatively, the *
c           control   can   be  viewed  as  attempting  to  provide  a *
c           calculated value of y at xend which is the exact  solution *
c           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
c           is proportional to tol.  (the norm  is  a  max  norm  with *
c           weights  that  depend on the error control strategy chosen *
c           by the user.  the default weight for the k-th component is *
c           1/max(1,abs(y(k))),  which therefore provides a mixture of *
c           absolute and relative error control.)                      *
c                                                                      *
c   ind  indicator - on initial entry ind must be set equal to  either *
c           1  or  2. if the user does not wish to use any options, he *
c           should set ind to 1 - all that remains for the user to  do *
c           then  is  to  declare c and w, and to specify nw. the user *
c           may also  select  various  options  on  initial  entry  by *
c           setting ind = 2 and initializing the first 9 components of *
c           c as described in the next section.  he may also  re-enter *
c           the  subroutine  with ind = 3 as mentioned again below. in *
c           any event, the subroutine returns with ind equal to        *
c              3 after a normal return                                 *
c              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
c              -1, -2, or -3 after an error condition (see below)      *
c                                                                      *
c     c  communications vector - the dimension must be greater than or *
c           equal to 24, unless option c(1) = 4 or 5 is used, in which *
c           case the dimension must be greater than or equal to n+30   *
c                                                                      *
c    nw  first dimension of workspace w -  must  be  greater  than  or *
c           equal to n                                                 *
c                                                                      *
c     w  workspace matrix - first dimension must be nw and second must *
c           be greater than or equal to 9                              *
c                                                                      *
c     the subroutine  will  normally  return  with  ind  =  3,  having *
c replaced the initial values of x and y with, respectively, the value *
c of xend and an approximation to y at xend.  the  subroutine  can  be *
c called  repeatedly  with new values of xend without having to change *
c any other argument.  however, changes in tol, or any of the  options *
c described below, may also be made on such a re-entry if desired.     *
c                                                                      *
c     three error returns are also possible, in which  case  x  and  y *
c will be the most recently accepted values -                          *
c     with ind = -3 the subroutine was unable  to  satisfy  the  error *
c        requirement  with a particular step-size that is less than or *
c        equal to hmin, which may mean that tol is too small           *
c     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
c        probably  means  that the requested tol (which is used in the *
c        calculation of hmin) is too small                             *
c     with ind = -1 the allowed maximum number of fcn evaluations  has *
c        been  exceeded,  but  this  can only occur if option c(7), as *
c        described in the next section, has been used                  *
c                                                                      *
c     there are several circumstances that will cause the calculations *
c to  be  terminated,  along with output of information that will help *
c the user determine the cause of  the  trouble.  these  circumstances *
c involve  entry with illegal or inconsistent values of the arguments, *
c such as attempting a normal  re-entry  without  first  changing  the *
c value of xend, or attempting to re-enter with ind less than zero.    *
c                                                                      *
c***********************************************************************
c                                                                      *
c     options - if the subroutine is entered with ind = 1, the first 9 *
c components of the communications vector are initialized to zero, and *
c the subroutine uses only default values  for  each  option.  if  the *
c subroutine  is  entered  with ind = 2, the user must specify each of *
c these 9 components - normally he would first set them all  to  zero, *
c and  then  make  non-zero  those  that  correspond to the particular *
c options he wishes to select. in any event, options may be changed on *
c re-entry  to  the  subroutine  -  but if the user changes any of the *
c options, or tol, in the course of a calculation he should be careful *
c about  how  such changes affect the subroutine - it may be better to *
c restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
c program  -  the information is available to the user, but should not *
c normally be changed by him.)                                         *
c                                                                      *
c  c(1)  error control indicator - the norm of the local error is  the *
c           max  norm  of  the  weighted  error  estimate  vector, the *
c           weights being determined according to the value of c(1) -  *
c              if c(1)=1 the weights are 1 (absolute error control)    *
c              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
c                 control)                                             *
c              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
c                 (relative  error  control,  unless abs(y(k)) is less *
c                 than the floor value, abs(c(2)) )                    *
c              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
c                 (here individual floor values are used)              *
c              if c(1)=5 the weights are 1/abs(c(k+30))                *
c              for all other values of c(1), including  c(1) = 0,  the *
c                 default  values  of  the  weights  are  taken  to be *
c                 1/max(1,abs(y(k))), as mentioned earlier             *
c           (in the two cases c(1) = 4 or 5 the user must declare  the *
c           dimension of c to be at least n+30 and must initialize the *
c           components c(31), c(32), ..., c(n+30).)                    *
c                                                                      *
c  c(2)  floor value - used when the indicator c(1) has the value 3    *
c                                                                      *
c  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
c           to be abs(c(3)) - otherwise it uses the default value      *
c              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
c           where dwarf is a very small positive  machine  number  and *
c           rreb is the relative roundoff error bound                  *
c                                                                      *
c  c(4)  hstart specification - if not zero, the subroutine  will  use *
c           an  initial  hmag equal to abs(c(4)), except of course for *
c           the restrictions imposed by hmin and hmax  -  otherwise it *
c           uses the default value of hmax*(tol)**(1/6)                *
c                                                                      *
c  c(5)  scale specification - this is intended to be a measure of the *
c           scale of the problem - larger values of scale tend to make *
c           the method more reliable, first  by  possibly  restricting *
c           hmax  (as  described  below) and second, by tightening the *
c           acceptance requirement - if c(5) is zero, a default  value *
c           of  1  is  used.  for  linear  homogeneous  problems  with *
c           constant coefficients, an appropriate value for scale is a *
c           norm  of  the  associated  matrix.  for other problems, an *
c           approximation to  an  average  value  of  a  norm  of  the *
c           jacobian along the trajectory may be appropriate           *
c                                                                      *
c  c(6)  hmax specification - four cases are possible                  *
c           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
c              min(abs(c(6)),2/abs(c(5)))                              *
c           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
c           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
c              2/abs(c(5))                                             *
c           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
c              of 2                                                    *
c                                                                      *
c  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
c           error  return with ind = -1 will be caused when the number *
c           of function evaluations exceeds abs(c(7))                  *
c                                                                      *
c  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
c           interrupt   the  calculations  after  it  has  chosen  its *
c           preliminary value of hmag, and just before choosing htrial *
c           and  xtrial  in  preparation for taking a step (htrial may *
c           differ from hmag in sign, and may  require  adjustment  if *
c           xend  is  near) - the subroutine returns with ind = 4, and *
c           will resume calculation at the point  of  interruption  if *
c           re-entered with ind = 4                                    *
c                                                                      *
c  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
c           interrupt   the  calculations  immediately  after  it  has *
c           decided whether or not to accept the result  of  the  most *
c           recent  trial step, with ind = 5 if it plans to accept, or *
c           ind = 6 if it plans to reject -  y(*)  is  the  previously *
c           accepted  result, while w(*,9) is the newly computed trial *
c           value, and w(*,2) is the unweighted error estimate vector. *
c           the  subroutine  will  resume calculations at the point of *
c           interruption on re-entry with ind = 5 or 6. (the user  may *
c           change ind in this case if he wishes, for example to force *
c           acceptance of a step that would otherwise be rejected,  or *
c           vice versa. he can also restart with ind = 1 or 2.)        *
c                                                                      *
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
c  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
c                                                                      *
c***********************************************************************
c                                                                      *
c  an overview of the program                                          *
c                                                                      *
c     begin initialization, parameter checking, interrupt re-entries   *
c  ......abort if ind out of range 1 to 6                              *
c  .     cases - initial entry, normal re-entry, interrupt re-entries  *
c  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
c  v........abort if n.gt.nw or tol.le.0                               *
c  .        if initial entry without options (ind .eq. 1)              *
c  .           set c(1) to c(9) equal to zero                          *
c  .        else initial entry with options (ind .eq. 2)               *
c  .           make c(1) to c(9) non-negative                          *
c  .           make floor values non-negative if they are to be used   *
c  .        end if                                                     *
c  .        initialize rreb, dwarf, prev xend, flag, counts            *
c  .     case 2 - normal re-entry (ind .eq. 3)                         *
c  .........abort if xend reached, and either x changed or xend not    *
c  .        re-initialize flag                                         *
c  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
c  v        transfer control to the appropriate re-entry point.......  *
c  .     end cases                                                  .  *
c  .  end initialization, etc.                                      .  *
c  .                                                                v  *
c  .  loop through the following 4 stages, once for each trial step .  *
c  .     stage 1 - prepare                                          .  *
c***********error return (with ind=-1) if no of fcn evals too great .  *
c  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
c  .        calc hmin, scale, hmax                                  .  *
c***********error return (with ind=-2) if hmin .gt. hmax            .  *
c  .        calc preliminary hmag                                   .  *
c***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
c  .        calc hmag, xtrial and htrial                            .  *
c  .     end stage 1                                                .  *
c  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
c  .     stage 3 - calc the error estimate                          .  *
c  .     stage 4 - make decisions                                   .  *
c  .        set ind=5 if step acceptable, else set ind=6            .  *
c***********interrupt no 2 if requested....................re-entry.v  *
c  .        if step accepted (ind .eq. 5)                              *
c  .           update x, y from xtrial, ytrial                         *
c  .           add 1 to no of successful steps                         *
c  .           set no of successive failures to zero                   *
c**************return(with ind=3, xend saved, flag set) if x .eq. xend *
c  .        else step not accepted (ind .eq. 6)                        *
c  .           add 1 to no of successive failures                      *
c**************error return (with ind=-3) if hmag .le. hmin            *
c  .        end if                                                     *
c  .     end stage 4                                                   *
c  .  end loop                                                         *
c  .                                                                   *
c  begin abort action                                                  *
c     output appropriate  message  about  stopping  the  calculations, *
c        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
c        previous xend,  no of  successful  steps,  no  of  successive *
c        failures, no of fcn evals, and the components of y            *
c     stop                                                             *
c  end abort action                                                    *
c                                                                      *
c***********************************************************************
*/

void dverk(dverkcomm *c, void (*fcn)(), double x, double *y, double xend, double tol, int *ind, int nw, void *fdata)
{
//     ******************************************************************
//     * begin initialization, parameter checking, interrupt re-entries *
//     ******************************************************************
//
//  ......abort if ind out of range 1 to 6

        int k;
        double temp;
//      int returnflag = 0;
        int n = c->n;
        double tout = x + c->ctrl->dt;
        
        double *w0 = (double *)calloc(n, sizeof(double));
        double *w1 = (double *)calloc(n, sizeof(double));
        double *w2 = (double *)calloc(n, sizeof(double));
        double *w3 = (double *)calloc(n, sizeof(double));
        double *w4 = (double *)calloc(n, sizeof(double));
        double *w5 = (double *)calloc(n, sizeof(double));
        double *w6 = (double *)calloc(n, sizeof(double));
        double *w7 = (double *)calloc(n, sizeof(double));
        double *w8 = (double *)calloc(n, sizeof(double));

        if ((*ind<1) || (*ind>6)) {
                goto L500;
        }
//        cases - initial entry, normal re-entry, interrupt re-entries
//         go to (5, 5, 45, 1111, 2222, 2222), ind
        switch (*ind) {
                case 1:
                        goto L5;
                        break;
                case 2:
                        goto L5;
                case 3:
                        goto L45;
                case 4:
                        goto L1111;
                case 5:
                        goto L2222;
                case 6:
                        goto L2222;
                default:
                        goto L500;
        }

//        case 1 - initial entry (ind .eq. 1 or 2)
//  .........abort if n.gt.nw or tol.le.0
L5:     if ((n>nw) || (tol<0.0)) {
                goto L500;
        }
        if (*ind==2) {
                goto L15;
        }
//              initial entry without options (ind .eq. 1)
//              set c(1) to c(9) equal to 0
/*      for(k = 0; k<9; k++) {
                c[k] = 0.0e0;
        } */
        c->errorcontrol = 0;
        c->floorvalue = 0.0;
        c->hminspec = 0.0;
        c->hstart = 0.0;
        c->scalespec = 0.0;
        c->hmaxspec = 0.0;
        c->maxevals = 0;
        c->interrupt1 = 0;
        c->interrupt2 = 0;
    goto L35;
L15:
//              initial entry with options (ind .eq. 2)
//              make c(1) to c(9) non-negative
//               do 20 k = 1, 9
/*      for(k=0; k<9; k++) {
        c[k] = fabs(c[k]);
        } */
        c->errorcontrol = abs(c->errorcontrol);
        c->floorvalue =fabs(c->floorvalue);
        c->hminspec = fabs(c->hminspec);
        c->hstart = fabs(c->hstart);
        c->scalespec = fabs(c->scalespec);
        c->hmaxspec = fabs(c->hmaxspec);
        c->maxevals = abs(c->maxevals);
        c->interrupt1 = abs(c->interrupt1);
        c->interrupt2 = abs(c->interrupt2);
L20:
//              make floor values non-negative if they are to be used
        if ((c->errorcontrol!=4) && (c->errorcontrol!=5)) {
                goto L30;
        }

//                  do 25 k = 1, n
        for(k=0; k<n; k++) {
                c->floorvalues[k] = fabs(c->floorvalues[k]);
        }
L30:
L35:
//           initialize rreb, dwarf, prev xend, flag, counts
        c->rreb = pow(2.0e0, -56);
        c->dwarf = 1.0e-35;
//           set previous xend initially to initial value of x
        c->prevxend = x;
//            do 40 k = 21, 24
/*      for(k=20; k<24; k++) {
                c[k] = 0.0e0;
        } */
        c->xendflag = 0;
        c->nsuccessful = 0;
        c->nfailed = 0;
        c->nevaluations = 0;

/*      if (c->ctrl->dt>0.0) {
                fcn(n, x, y, w0, fdata);
                c->nevaluations = c->nevaluations + 1;
                c->ctrl->output(n, y, x);
                tout += c->ctrl->dt;
        } */
                
        goto L50;
//        case 2 - normal re-entry (ind .eq. 3)
//  .........abort if xend reached, and either x changed or xend not
L45:
        if ((c->xendflag!=0) && (x!=c->prevxend || xend==c->prevxend)) {
                goto L500;
        }
//           re-initialize flag
        c->xendflag = 0;
        goto L50;
//        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
//           transfer control to the appropriate re-entry point..........
//           this has already been handled by the computed go to        .
//        end cases                                                     v
L50:
/*    end initialization, etc.
c
c     ******************************************************************
c     * loop through the following 4 stages, once for each trial  step *
c     * until the occurrence of one of the following                   *
c     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
c     *        stage 4                                                 *
c     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
c     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
c     *        requested, in stage 1 or stage 4                        *
c     ****************************************************************** */

L99999:
        if ((x>tout) && (c->ctrl->dt>0.0)) {
                c->ctrl->output(c->ctrl, n, y, x);
                tout += c->ctrl->dt;
        }

/*        ***************************************************************
c        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
c        * and some parameter  checking,  and  end  up  with  suitable *
c        * values of hmag, xtrial and htrial in preparation for taking *
c        * an integration step.                                        *
c        *************************************************************** */

//      error return (with ind=-1) if no of fcn evals too great
        if (c->maxevals==0 || c->nevaluations<c->maxevals) {
                goto L100;
        }
        *ind = -1;
    goto L5000;
        
L100:

//      calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (*ind==6) {
                goto L105;
        }
        fcn(n, x, y, w0, fdata);
    c->nevaluations = c->nevaluations + 1;
//      c->ctrl->output(n, y, x);
//      tout += c->ctrl->dt;
L105:

//      calculate hmin - use default unless value prescribed
        c->hmin = c->hminspec;
        if (c->hminspec != 0.0e0) goto L165;
//              calculate default value of hmin
//              first calculate weighted norm y - c(12) - as specified
//              by the error control indicator c(1)
        temp = 0.0e0;
        if (c->errorcontrol != 1) goto L115;
//                 absolute error control - weights are 1
//                  do 110 k = 1, n
        for(k=0; k<n; k++) {
                temp = (temp>fabs(y[k]))?temp:fabs(y[k]);
        }
        c->weightednormy = temp;
        goto L160;
L115:          
        if (c->errorcontrol != 2) goto L120;
//                 relative error control - weights are 1/dabs(y(k)) so
//                 weighted norm y is 1
    c->weightednormy = 1.0e0;
    goto L160;
L120:          
        if (c->errorcontrol != 3) goto L130;
//                 weights are 1/max(c(2),abs(y(k)))
//                  do 125 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(y[k])/c->floorvalue);
        }
        c->weightednormy = dmin1(temp, 1.0e0);
        goto L160;
L130:          
        if (c->errorcontrol != 4) goto L140;
//                 weights are 1/max(c(k+30),abs(y(k)))
//                  do 135 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(y[k])/c->floorvalues[k]);
        }
        c->weightednormy = dmin1(temp, 1.0e0);
        goto L160;
L140:          
        if (c->errorcontrol != 5) goto L150;
//                 weights are 1/c(k+30)
//                  do 145 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(y[k])/c->floorvalues[k]);
        }
        c->weightednormy = temp;
        goto L160;
L150:
//                 default case - weights are 1/max(1,abs(y(k)))
//                  do 155 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(y[k]));
        }
        c->weightednormy = dmin1(temp, 1.0e0);
L160:
        c->hmin = 10.0e0*dmax1(c->dwarf, c->rreb*dmax1(c->weightednormy/tol,fabs(x)));
L165:

//           calculate scale - use default unless value prescribed
        c->scale = c->scalespec;
        if (c->scalespec == 0.0e0) c->scale = 1.0e0;

//           calculate hmax - consider 4 cases
//          case 1 both hmax and scale prescribed
        if (c->hmaxspec!=0.0e0 && c->scalespec!=0.0e0) c->hmax = dmin1(c->hmaxspec, 2.0e0/c->scalespec);
//          case 2 - hmax prescribed, but scale not
        if (c->hmaxspec!=0.0e0 && c->scalespec==0.0e0) c->hmax = c->hmaxspec;
//           case 3 - hmax not prescribed, but scale is
        if (c->hmaxspec==0.0e0 && c->scalespec!=0.0e0) c->hmax = 2.0e0/c->scalespec;
//           case 4 - neither hmax nor scale is provided
        if (c->hmaxspec==0.0e0 && c->scalespec==0.0e0) c->hmax = 2.0e0;

//***********error return (with ind=-2) if hmin .gt. hmax
        if (c->hmin < c->hmax) goto L170;
        *ind = -2;
        goto L5000;
L170:

//           calculate preliminary hmag - consider 3 cases
        if (*ind > 2) goto L175;
//           case 1 - initial entry - use prescribed value of hstart, if
//              any, else default
        c->hmag = c->hstart;
        if (c->hstart == 0.0e0) c->hmag = c->hmax*pow(tol, (1.0/6.0));
        goto L185;
L175:
        if (c->nfailed > 1) goto L180;
//           case 2 - after a successful step, or at most  one  failure,
//              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
//              overflow. then avoid reduction by more than half.
        temp = 2.0e0*c->hmag;
        if (tol < pow(2.0e0/0.9e0, 6)*c->est) temp = 0.9e0*tol/pow(c->est, (1.0/6.0))*c->hmag;
        c->hmag = dmax1(temp, 0.5e0*c->hmag);
        goto L185;
L180:
//           case 3 - after two or more successive failures
        c->hmag = 0.5e0*c->hmag;
L185:

//           check against hmax
        c->hmag = dmin1(c->hmag, c->hmax);

//           check against hmin
        c->hmag = dmax1(c->hmag, c->hmin);

//***********interrupt no 1 (with ind=4) if requested
        if (c->interrupt1 == 0) goto L1111;
        *ind = 4;
        goto L5000;
//           resume here on re-entry with ind .eq. 4   ........re-entry..
L1111:

//           calculate hmag, xtrial - depending on preliminary hmag, xend
        if (c->hmag >= fabs(xend - x)) goto L190;
//              do not step more than half way to xend
        c->hmag = dmin1(c->hmag, 0.5e0*fabs(xend - x));
        c->xtrial = x + dsign(c->hmag, xend - x);
    goto L195;
L190:
//              hit xend exactly
        c->hmag = fabs(xend - x);
    c->xtrial = xend;
L195:

//           calculate htrial
        c->htrial = c->xtrial - x;
/*
c        end stage 1
c
c        ***************************************************************
c        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
c        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
c        * stage 3. w(*,9) is temporary storage until finally it holds *
c        * ytrial.                                                     *
c        *************************************************************** */
//      printf("%15.6e", c[16]);
/*      printf("%8.5f", x);
        for(k=0; k<n; k++) {
                printf(" %15.6e", y[k]);
        }
        printf("\n");  */

//      Here!
/*      if ((x>tout) && (c->ctrl->dt>0.0)) {
                c->ctrl->output(c->ctrl, n, y, x);
                tout += c->ctrl->dt;
        } */

        temp = c->htrial/1398169080000.0e0;

//            do 200 k = 1, n
        for(k=0; k<n; k++) {
                w8[k] = y[k] + temp*w0[k]*233028180000.0e0;
        }
    fcn(n, x + c->htrial/6.0e0, w8, w0, fdata);

//            do 205 k = 1, n
        for(k=0; k<n; k++) {
                w8[k] = y[k] + temp*(w0[k]*74569017600.0e0
                                           + w1[k]*298276070400.0e0  );
        }
    fcn(n, x + c->htrial*(4.0e0/15.0e0), w8, w2, fdata);

//            do 210 k = 1, n
        for(k=0; k<n; k++) {
                w8[k] = y[k] + temp*(w0[k]*1165140900000.0e0
                                           - w1[k]*3728450880000.0e0
                                           + w2[k]*3495422700000.0e0 );
        }
    fcn(n, x + c->htrial*(2.0e0/3.0e0), w8, w3, fdata);

//            do 215 k = 1, n
        for(k=0; k<n; k++) {
                w8[k] = y[k] + temp*(-w0[k]*3604654659375.0e0
                                                          + w1[k]*12816549900000.0e0 
                                                          - w2[k]*9284716546875.0e0 
                                                          + w3[k]*1237962206250.0e0 );
        }
    fcn(n, x + c->htrial*(5.0e0/6.0e0), w8, w4, fdata);

//            do 220 k = 1, n
        for(k=0; k<n; k++) {
                w8[k] = y[k] + temp*(w0[k]*3355605792000.0e0 
                                     - w1[k]*11185352640000.0e0 
                                                         + w2[k]*9172628850000.0e0 
                                                         - w3[k]*427218330000.0e0 
                                                         + w4[k]*482505408000.0e0);
        }
    fcn(n, x + c->htrial, w8, w5, fdata);

//            do 225 k = 1, n
        for(k=0; k<n; k++) {
                w8[k] = y[k] + temp*(-w0[k]*770204740536.0e0
                                     + w1[k]*2311639545600.0e0
                                     - w2[k]*1322092233000.0e0
                                     - w3[k]*453006781920.0e0
                                     + w4[k]*326875481856.0e0);
        }
        fcn(n, x + c->htrial/15.0e0, w8, w6, fdata);

//            do 230 k = 1, n
        for(k=0; k<n; k++) {
                w8[k] = y[k] + temp*(w0[k]*2845924389000.0e0
                                    - w1[k]*9754668000000.0e0
                                    + w2[k]*7897110375000.0e0
                                    - w3[k]*192082660000.0e0
                                    + w4[k]*400298976000.0e0
                                    + w6[k]*201586000000.0e0  );
        }
    fcn(n, x + c->htrial, w8, w7, fdata);

//           calculate ytrial, the extrapolated approximation and store
//              in w(*,9)
//            do 235 k = 1, n
        for(k=0; k<n; k++) {
                w8[k] = y[k] + temp*(w0[k]*104862681000.0e0
                                     + w2[k]*545186250000.0e0
                                     + w3[k]*446637345000.0e0
                                     + w4[k]*188806464000.0e0
                                     + w6[k]*15076875000.0e0
                                     + w7[k]*97599465000.0e0   );
        }

//           add 7 to the no of fcn evals
        c->nevaluations = c->nevaluations + 7;
/*
c        end stage 2
c
c        ***************************************************************
c        * stage 3 - calculate the error estimate est. first calculate *
c        * the  unweighted  absolute  error  estimate vector (per unit *
c        * step) for the unextrapolated approximation and store it  in *
c        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
c        * specified by the error  control  indicator  c(1).  finally, *
c        * modify  this result to produce est, the error estimate (per *
c        * unit step) for the extrapolated approximation ytrial.       *
c        ***************************************************************
c
c           calculate the unweighted absolute error estimate vector */
//            do 300 k = 1, n
        for(k=0; k<n; k++) {
                w1[k] = (w0[k]*8738556750.0e0
                         + w2[k]*9735468750.0e0
                         - w3[k]*9709507500.0e0
                         + w4[k]*8582112000.0e0
                         + w5[k]*95329710000.0e0
                         - w6[k]*15076875000.0e0
                         - w7[k]*97599465000.0e0)/1398169080000.0e0;
        }

//           calculate the weighted max norm of w(*,2) as specified by
//           the error control indicator c(1)
        temp = 0.0e0;
    if (c->errorcontrol != 1) goto L310;
//              absolute error control
//               do 305 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp,fabs(w1[k]));
        }
    goto L360;
L310:       
        if (c->errorcontrol != 2) goto L320;
//              relative error control
//               do 315 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(w1[k]/y[k]));
        }
    goto L360;
L320:
        if (c->errorcontrol != 3) goto L330;
//              weights are 1/max(c(2),abs(y(k)))
//               do 325 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(w1[k]) / dmax1(c->floorvalue, fabs(y[k])) );
        }
    goto L360;
L330:
        if (c->errorcontrol != 4) goto L340;
//              weights are 1/max(c(k+30),abs(y(k)))
//               do 335 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(w1[k]) / dmax1(c->floorvalues[k], fabs(y[k])) );
        }
    goto L360;
L340:
    if (c->errorcontrol != 5) goto L350;
//              weights are 1/c(k+30)
//               do 345 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(w1[k]/c->floorvalues[k]));
        }
    goto L360;
L350:
//              default case - weights are 1/max(1,abs(y(k)))
//               do 355 k = 1, n
        for(k=0; k<n; k++) {
                temp = dmax1(temp, fabs(w1[k]) / dmax1(1.0e0, fabs(y[k])) );
        }
L360:

/*           calculate est - (the weighted max norm of w(*,2))*hmag*scale
c              - est is intended to be a measure of the error  per  unit
c              step in ytrial */
        c->est = temp*c->hmag*c->scale;
/*
c        end stage 3
c
c        ***************************************************************
c        * stage 4 - make decisions.                                   *
c        ***************************************************************
c
c           set ind=5 if step acceptable, else set ind=6 */
        *ind = 5;
        if (c->est > tol) *ind = 6;

//**********interrupt no 2 if requested
        if (c->interrupt2 == 0) goto L2222;
        goto L5000;
//          resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
L2222:

        if (*ind == 6) goto L410;
//              step accepted (ind .eq. 5), so update x, y from xtrial,
//                 ytrial, add 1 to the no of successful steps, and set
//                 the no of successive failures to zero
        x = c->xtrial;
//               do 400 k = 1, n
        for(k=0; k<n; k++) {
                y[k] = w8[k];
        }
        (c->nsuccessful)++;
    c->nfailed = 0;
//**************return(with ind=3, xend saved, flag set) if x .eq. xend
        if (x != xend) goto L405;
    *ind = 3;
    c->prevxend = xend;
    c->xendflag = 1;
    goto L5000;
L405:
        goto L420;
L410:
//              step not accepted (ind .eq. 6), so add 1 to the no of
//                 successive failures
        (c->nfailed)++;// = c[22] + 1.0e0;
//**************error return (with ind=-3) if hmag .le. hmin
        if (c->hmag > c->hmin) goto L415;
    *ind = -3;
    goto L5000;
L415:
L420:

//        end stage 4

        goto L99999;
//     end loop

//  begin abort action
L500:

//      write(6,505) ind, tol, x, n, c(13), xend, nw, c(16), c(20),
//     +      c(22), c(23), c(24), (y(k), k = 1, n)
//
        printf("Computation stopped on dverk: %4d %15.6e %15.6e %4d %15.6e %15.6e\n", *ind,tol,x,n,c->weightednormy,xend);

/*  505 format( /// 1h0, 58hcomputation stopped in dverk with the followin
c     +g values -
c     +   / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =,
c     +          1pd22.15
c     +   / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =,
c     +          1pd22.15
c     +   / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =,
c     +          1pd22.15
c     +   / 1h0, 14x, 27hno of successful steps    =, 0pf8.0
c     +   / 1h , 14x, 27hno of successive failures =, 0pf8.0
c     +   / 1h , 14x, 27hno of function evals      =, 0pf8.0
c     +   / 1h0, 23hthe components of y are
c     +   // (1h , 1p5d24.15)                                           )*/

//505   format(///'Computation stopped on dverk:'/
//     +                  i4,2d15.6,i4,2d14.5///)
L5000:
        if (c->ctrl->dt>0.0) {
                fcn(n, x, y, w0, fdata);
                c->nevaluations = c->nevaluations + 1;
                c->ctrl->output(c->ctrl, n, y, x);
//              tout += c->ctrl->dt;
        } 
        
        free(w0);
        free(w1);
        free(w2);
        free(w3);
        free(w4);
        free(w5);
        free(w6);
        free(w7);
        free(w8);
        
    return;

}