#ifndef FILTER_H_
#define FILTER_H_
//#define PI 3.14159
#include "util.h"
#define px(x,y,lengthx, lengthy) (x+y*lengthx)

void createFilter(int radius, double sigma, double** &gKernel);

void conv(double* in, int rows, int cols, double **kernel,int kRows, int kCols, double** &out);

double* conv2(double* in, int rows, int cols, double **kernel,int kRows, int kCols);
#endif