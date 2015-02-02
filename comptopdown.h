//#complearn.h
#ifndef COMPTOPDOWN_H_
#define COMPTOPDOWN_H_
#include "util.h"
#include <assert.h>
/*
Data* comptopdown(double *morients, int imh, int imw, struct Parts *model, int numlayer, int numsub, 
        int numOrient, int* radius, double discount, double *param);*/

/* ----------------------------------------------------------
 *                  main function 
 * ----------------------------------------------------------
 * /

/*single image inference*/
Data* comptopdown(double *morients,int imh, int imw, struct Parts *model, int numlayer,int numsub, int numOrient,double discount, double* param);

/*single image with model rotation*/
void comptopdown(Data *output, unsigned char* morients, int imh, int imw, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, double rotation);
/*process each image*/
void comptopdown_fast(Data *output, unsigned char* morients, int imh, int imw, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, double rotation);

void comptopdown_full(Data *output, unsigned char* morients, int imh, int imw, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, double rotation);
/*batch processing*/
Data* comptopdown_ex(struct Data *examples, int ndata, struct Parts *model, int numlayer, int numsub, int numOrient, double discount, double *radjust);
/*batch processing with rotation*/
Data* comptopdown_ex(struct Data *examples, int ndata, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, int rotation);
//fast one
Data* comptopdown_ex2(struct Data *examples, int ndata, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, int rotation);

/* -----------------------------------------------------
 *              matching subparts
 * -----------------------------------------------------
 */
bool topdownmatch(double *pinstances, int ninst, int col, int o1, int r, int c, int *dyx, int *o2list, int numdyx, Parts *model, int numlayer, 
        int numOrient, int numsub, double *radjust, double discount, Node *point);


/*matching subparts with rotation on the fly*/
bool topdownmatch(double *pinstances, int ninst, int col, int o1, int r, int c, int *dyx, int *o2list, int numdyx, Parts *model, int numlayer, 
        int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point);
/*matching with rotation old version*/
bool topdownmatch_old(double *pinstances, int ninst, int col, int o1, int r, int c, int *dyx, int *o2list,int numdyx,
        Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point);

/*matching only children*/
bool topdownmatch_fast(double *pinstances, int ninst, int pDim, int *pdist, int iinst, struct array1d idx, bool *visited, Parts *model, 
        int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point);

/*matching both anchor and children*/
bool topdownmatch_full(double *pinstances, int ninst, int col, int *pdist, int iinst, struct array1d idx, bool *visited,
        Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point);


/* ----------------------------------------------------------
 *                        utility function
 * ----------------------------------------------------------
 * /

/*get the number of childs from top-down*/
int getnumchild(int centerid, struct Parts *model, int numlayer, int numsub);
//struct array1d matchmid(int o1, struct Parts *model,int numlayer, int numsub);

int id2orient(int centerid, struct Parts *model, int numlayer, int numsub);

double* drawsubparts(int *img, int imh, int imw, int y, int x,   struct Parts *model,int centerid, int numsub, int numOrient, int numlayer, int **allFilter, int h, int w);

struct pTypes* mappingParts(struct Parts *model, int numOrient, int numlayer, int numsub, struct array1d *tdlinks);
/*mapping parts with rotation on the fly*/
struct pTypes* mappingParts(struct Parts *model, int numOrient, int numlayer, int numsub, int rotation, struct array1d *tdlinks);
/*flatten all layer parts with orientations*/
struct fParts* flattenparts(struct Parts *model, int numOrient, int numlayer, int numsub, int rotation);

struct array1d unique(struct array1d d, int dim);

void initialize(int *flag, int n);
void initialize(bool *visited, int n);

int vmax(double *vector, int n, double *score);

struct Node cat(struct Node pchild, struct Node cchild);

struct Node cat(struct Node pchild, struct Node cchild, int pDim);

void cat(struct Node *pchild, struct Node *cchild, int pDim);

/* ---------------------------------------------------------------------
 *                          exact matching here
 * ---------------------------------------------------------------------
 * /

/*TOP DOWN MATCHING WITH ROTATION, deleting matched point in the next recursive step*/
bool exactmatch(double *pinstances, int ninst, int col, int *pdist, int iinst, struct array1d idx, bool *visited,
        Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point);

/*match specific subchildren from top down hierarchy*/
bool hiertpdmatch(double *pinstances, int ninst, int col, int *pdist, int iinst, struct array1d idx, bool *visited,
        Parts *model, int pid, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point);

#endif