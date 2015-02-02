//#complearn.h
#ifndef COMPINFERENCE_H_
#define COMPINFERENCE_H_
#include "util.h"
#include "checkparts.h"


Data* compinference(struct Data *examples, int numdata, struct Parts *model,int numlayer, int numsub,  int numOrient, int radius, 
        double percent, double *param, double discount);
/*
Data* compinference(struct Data *examples, int numdata, struct Parts *model,int numlayer, int numsub,  int numOrient, int radius, 
        double percent, double *param, double discount, bool flag);
*/

Data* compinference(struct Data *examples, int ndata, struct Parts *model,int numlayer, int numscale, int numsub,  int numOrient, int radius, 
        double percent, double *param, double discount);


struct double_array1d matchingparts(int o1, int r, int c, int *dyx, int *o2list,int numdyx,struct Parts * model,
        int numlayer,int subindex, double discount);

struct array1d matchmid(int o1, struct Parts *model,int numlayer, int numsub);

#endif