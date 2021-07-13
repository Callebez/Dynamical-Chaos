#ifndef PLOTTING_H
#define PLOTTING_H
#include <string>
#include<iostream>
void plot2D(const char* fileName,const char* inputFile, const char* titles, const char* optionalFlagsGnuPlot);
void plot3D(const char* fileName,const char* inputFile, const char* titles, const char* optionalFlagsGnuPlot);
void plotAnimate2D(const char* fileName,const char* inputFile, int delay);
void plotAnimate3D(const char* fileName,const char* inputFile, int delay);
#endif