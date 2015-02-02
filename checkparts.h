//#complearn.h
#ifndef CHECKPARTS_H_
#define CHECKPARTS_H_
#include "util.h"
/*matching parts with orientation and relative positions*/
struct array1d checksubparts(int o1, int r, int c, int* dyx, int* o2list, int numdyx, double* ptypes, int* nptypes, 
        int pDim,int numlayer, int subindex, double discount, struct array1d *lutlinks);
//find parts whoes anchorid equals to o1
struct array1d searchanchorid(int o1, double* ptypes, int* nptypes, int pDim);
/*
 *match from the model
struct array1d checksubparts(int o1, int r, int c, int* dyx, int* o2list, int numdyx, struct pTypes *ptypes, int* nptypes, 
        int numlayer, int subindex, double discount);*/
struct array1d checksubparts(int o1, int r, int c, int* dyx, int* o2list, int numdyx, struct pTypes *ptypes, int* nptypes, 
        int numlayer, int subindex, double discount, struct array1d *lutlinks);
struct array1d searchanchorid(int o1, struct pTypes *ptypes, int *nptypes, int subindex);

double* selection(int* dyx, int numdyx, int dyxDim, int* indicator);
/*suppression straight line by setting weights*/
double * setweights(struct pTypes *ptypes, int *nptypes, int subindex, double rate);
double * setweights(double *ptypes, int *nptypes,  int pDim, int subindex, double rate);

// adjust weight for two subcomposition
double adjweight(struct pTypes *ptypes, int *nptypes, int subindex, int anchorid, int cid, double rate);
double adjweight(double *ptypes, int *nptypes, int pDim, int subindex, int anchorid, int cid, double rate);

#endif
