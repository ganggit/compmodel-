//#complearn.h
#ifndef COMPLEARN_H_
#define COMPLEARN_H_
#include "util.h"

void learningcomposition(struct Data *examples, int numdata, struct Parts *model,int numlayer, int numsub,  int numOrient, int radius, 
        double percent, double *param, double discount, bool flag);

void initmask(int *img, int imh, int imw);
void initmask3(int *img, int imh, int imw, int pDim);
void buildmask(int *img, int imh, int imw, double *pinstances, int ninst, int col);
void buildmask3(int *img, int imh, int imw, int pDim, double *pinstances, int ninst, int col);
struct array1d lookneighs(int y, int x, int* img, int imh, int imw, int radius, int pre_radius);
struct array1d lookneighs3(int y, int x, int *img, int imh, int imw, int pDim, int radius, int pre_radius);

/*get the points in the histogram*/
int* getneighs(int y, int x, double *img, int imh, int imw, int radius, int *numpoint);
#endif