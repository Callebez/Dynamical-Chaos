#ifndef PLOTTING_H
#define PLOTTING_H
#include <string>
#include<iostream>
void plot2D(const char* fileName,const char* inputFile,const char* systemName, const char* keys);
void plot3D(const char* fileName,const char* inputFile,const char* systemName, const char* keys);
void plotAnimate2D(const char* fileName,const char* inputFile, int delay);
void plotAnimate3D(const char* fileName,const char* inputFile, int delay);
void plotBifucation(const char* fileName, const char* systemName, const char* paramName, const char* optionalFlagsGnuPlot);
void plotBifucationAndLyapunovExp(const char* biffFileName,const char* LyapunovFileName, const char* systemName, const char* paramName, const char* optionalFlagsGnuPlot);

#endif