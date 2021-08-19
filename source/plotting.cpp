#include "../include/plotting.hpp"

void plot2D(const char* fileName,const char* inputFile,const char* systemName,const char* keys)
{
    FILE *gnupipe = popen("gnuplot -persist", "w");
    if(gnupipe)
    {
        fprintf(gnupipe,"set terminal pngcairo enhanced size 1080,720 font \'Heveltica, 15\' \n");
        fprintf(gnupipe,"set title \"%s \" font \'Helvetica, 16\' \n ",systemName);
        fprintf(gnupipe,"set tics font \'Helvetica, 14\' \n" );
        fprintf(gnupipe,"set border lw 3 \n " );
        fprintf(gnupipe,"set key opaque \n" );
        fprintf(gnupipe,"set output \'./outputs/images/%s.png\' \n", fileName);
        fprintf(gnupipe,"plot \'./outputs/txt/%s.dat\' w lines lw 0.5 lc 7 title \'%s\'", fileName ,keys);
    }
}
void plot3D(const char* fileName,const char* inputFile,const char* systemName, const char* keys)
{
    FILE *gnupipe = popen("gnuplot -persist", "w");
    if(gnupipe)
    {
        fprintf(gnupipe,"set terminal pngcairo enhanced size 1080,720 font \'Heveltica, 15\' \n");
        fprintf(gnupipe,"set title \"%s \" font \'Helvetica, 16\' \n ", systemName);
        fprintf(gnupipe,"set tics font \'Helvetica, 14\' \n" );
        fprintf(gnupipe,"set border lw 3 \n " );
        fprintf(gnupipe,"set key opaque \n" );
        fprintf(gnupipe,"set output \'./outputs/images/%s.png\' \n", fileName);
        fprintf(gnupipe,"splot \'./outputs/txt/%s.dat\' w lines lw 0.5 lc 7 title \'%s\'", fileName ,keys);
    }

}
void plotAnimate3D(const char* fileName,const char* inputFile, int delay)
{
    FILE *gnupipe = popen("gnuplot -persist", "w");
    if(gnupipe)
    {
        fprintf(gnupipe,"set terminal gif animate delay %i\n", delay);
        fprintf(gnupipe,"stats \'%s\' name  \'A\' \n", inputFile);
        fprintf(gnupipe,"set output \'%s\'\n", fileName);
        fprintf(gnupipe,"do for [i=0:int(A_blocks-1)]{ splot \"%s\" index i w lines lw 2 lc 7 }\n", inputFile);
        fprintf(gnupipe,"\n");
    }
}
void plotAnimate2D(const char* fileName,const char* inputFile, int delay)
{
    FILE *gnupipe = popen("gnuplot -persist", "w");
    if(gnupipe)
    {
        fprintf(gnupipe,"set terminal gif animate delay %i\n", delay);
        fprintf(gnupipe,"stats \'%s\' name  \'A\' \n", inputFile);
        fprintf(gnupipe,"set output \'%s\'\n", fileName);
        fprintf(gnupipe, "set xrange[0:20]\n set yrange[-10:10]\n");
        fprintf(gnupipe,"do for [i=0:int(A_blocks-1)]{ plot \"%s\" index i w lines lw 2 lc 7 }\n", inputFile);
        fprintf(gnupipe,"\n");
    }
}
void plotBifucation(const char* fileName, const char* systemName, const char* paramName, const char* optionalFlagsGnuPlot)
{

    FILE *gnupipe = popen("gnuplot -persist", "w");
    if(gnupipe)
    {
        fprintf(gnupipe,"set terminal pngcairo enhanced size 1080,720 font \'Heveltica, 15\' \n");
        fprintf(gnupipe,"set title \"Bifucation diagram in x  \\n(%s) \" font \'Helvetica, 16\' \n ",systemName);
        fprintf(gnupipe,"set xlabel \'%s\' enhanced font \'Helvetica, 14\' \n", paramName);
        fprintf(gnupipe,"set ylabel \' X\' font \"Helvetica, 14\"\n");
        fprintf(gnupipe,"set tics font \'Helvetica, 14\' \n" );
        fprintf(gnupipe,"set border lw 3 \n " );
        fprintf(gnupipe,"set key opaque \n" );
        // fprintf(gnupipe,"%s \n", optionalFlagsGnuPlot);
        fprintf(gnupipe,"set output \'./outputs/images/%s.png\' \n", fileName);
        fprintf(gnupipe,"plot \'./outputs/txt/%smax.dat\' using 1:2  pt 0 ps 0.75 title \'max points\',  \'./outputs/txt/%smin.dat\'  using 1:2  pt 0 ps 0.75 title \'min points\'  \n", fileName ,fileName);
    }
}
void plotBifucationAndLyapunovExp(const char* biffFileName,const char* LyapunovFileName, const char* systemName, const char* paramName, const char* optionalFlagsGnuPlot)
{

    FILE *gnupipe = popen("gnuplot -persist", "w");
    if(gnupipe)
    {
        fprintf(gnupipe,"set terminal pngcairo enhanced size 1080,720 font \'Heveltica, 15\' \n");
        fprintf(gnupipe,"set title \"Bifucation diagram in x and Lyapunov Exponents  \\n(%s) \" font \'Helvetica, 16\' \n ",systemName);
        fprintf(gnupipe,"set xlabel \'%s\' enhanced font \'Helvetica, 14\' \n", paramName);
        fprintf(gnupipe,"set ylabel \' X\' font \"Helvetica, 14\"\n");
        fprintf(gnupipe,"set tics font \'Helvetica, 14\' \n" );
        fprintf(gnupipe,"set border lw 3 \n " );
        fprintf(gnupipe,"set key opaque \n" );
        fprintf(gnupipe,"set grid \n" );
        fprintf(gnupipe,"set yrange[-20:20] \n" );

        // fprintf(gnupipe,"%s \n", optionalFlagsGnuPlot);
        fprintf(gnupipe,"set output \'./outputs/images/%s.png\' \n", biffFileName);
        fprintf(gnupipe,"plot \'./outputs/txt/%smax.dat\' using 1:2  pt 0 ps 0.75 title \'max points\', \
               \'./outputs/txt/%smin.dat\'  using 1:2  pt 0 ps 0.75 title \'min points\', \
               \'./outputs/txt/%s1.dat\'  using 1:2 w lines lw 2 title \'Lyapunov exp\' ,\
               \'./outputs/txt/%s2.dat\'  using 1:2 w lines lw 2 title \'Lyapunov exp\' ,\
               \'./outputs/txt/%s3.dat\'  using 1:2 w lines lw 2 title \'Lyapunov exp\' \
                 \n", biffFileName ,biffFileName, LyapunovFileName,LyapunovFileName, LyapunovFileName);
    }
}