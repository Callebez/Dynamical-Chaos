#include "../include/plotting.hpp"

void plot2D(const char* fileName,const char* inputFile, const char* titles, const char* optionalFlagsGnuPlot)
{
    FILE *gnupipe = popen("gnuplot -persist", "w");
    if(gnupipe)
    {
        
        fprintf(gnupipe,"%s\n", optionalFlagsGnuPlot);
        fprintf(gnupipe,"set output \'%s.jpeg\'\n", fileName);
        fprintf(gnupipe,"plot \'%s\'", inputFile, "lw 3 lc 7 title \'%s\'\n", titles);
    }
}
void plot3D(const char* fileName,const char* inputFile, const char* titles, const char* optionalFlagsGnuPlot)
{
    FILE *gnupipe = fopen("gnuplot -persist", "w");
    if(gnupipe)
    {
        
        fprintf(gnupipe,"%s", optionalFlagsGnuPlot);
        fprintf(gnupipe,"set output \'%s.jpeg\'", fileName);
        fprintf(gnupipe,"splot \'%s\'", inputFile, "lw 2 lc 7 title \'%s\'", inputFile, titles);
    }
}
void plotAnimate(const char* fileName,const char* inputFile, int delay)
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