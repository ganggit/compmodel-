#ifndef UTIL_H_
#define UTIL_H_

#define MAX_LEN 5000
#define PI 3.14159265
#define INF 99999999
#define MAX_NUM_PARTS 2000
#define EXT_NUM_PARTS 10000
#define MAX_POINTS 100
#define MDIM 2
//#define MIN(a,b) = ((a<b)? a:b)
//#define MAX(a,b) = (a>b? a:b)
//#define numsub 2
#include <mex.h>
#include <string>
#include <math.h>
#include <time.h>
//#include <mcheck.h>
#include "filter.h"

/*address*/
#define px2(x, y, lengthx, lengthy) ((x) + (y)*(lengthx) - 0) 
#define px3(x, y, z, lengthx, lengthy, lengthz) ((x) + ((y)*(lengthx)) + ((z)*(lengthx)*(lengthy)) - 0) 
//#define px3(x, y, z, lengthx, lengthy, lengthz) ((z) + ((y)*(lengthz)) + ((x)*(lengthy)*(lengthz))) 
/*easy function
#define ABS(x) ((x)>0? (x):(-(x)))*/
#define MAX(x, y) ((x)>(y)? (x):(y))
#define MIN(x, y) ((x)<(y)? (x):(y))


/* define a struct to store matching index */

struct double_array1d 
{
   double a[MAX_LEN];
   double length;
};


struct array1d 
{
   int a[MAX_LEN];
   int length;
};

struct cellData
{
    double *a;
    int length;
};

//flatten parts model
struct fParts{
    double_array1d dyx[MAX_NUM_PARTS+EXT_NUM_PARTS]; // parts' ID and relative postions
    array1d *tdlinks; //bottom up links for fast matching
    int ntypes; // the number of parts
};


//for multiple subparts
struct pTypes{
 double_array1d dyx[MAX_NUM_PARTS+EXT_NUM_PARTS]; // parts' ID and relative postions
 double_array1d p; //probability for the parts
 int ntypes;
};

struct Parts{
    pTypes *ptypes_; // part types
    int *ntypes_; //fast pointer to the number of parts
    array1d *plinks_; //bottom up links for fast matching
    int radius;
    
};

/*row col of elements */
struct sPinstances{
    int row;
    int col;
    double *elements;
        
};
/*cell and num of cell*/
struct sLambda{
    
    //array1d elements[MAX_NUM_PARTS];
    cellData *elements;
    int numel;
    
};

struct Data{
    
    int imh;
    int imw;
    //std::string imname;
    char *imname;
    double *mask;
    double *morients;
    sPinstances *pinstances;
    sLambda *lambda;
    int numlayer;
};


struct Node{
double score;
int compid;
int a[4*MAX_LEN];
int length;
//struct Node *next;
};


/* Generating integer vector */
int *int_vector(int n);
/* Generating integer matrix */
int **int_matrix(int m, int n);

/* Generating doulbe vector */
double *double_vector(int n);
/* Generating double matrix */
double **double_matrix(int m, int n);


/* Free matrix space */
void free_matrix(void **mat, int m, int n);

/*dilate images*/
//int** dilate(int** image);
void imdilate(int *image, int imw, int imh);
/*check any intersection between image*/
double imoverlap(int *im1, int *im2, int imh, int imw);


void swap(double **a, int dim, int left, int right);
/*quick sort for matrix*/
int partition(double **a,int dim, int low, int high);
void quicksort(double **a,int dim, int low, int high);

void quickSort2(double arr[], int left, int right);

//bubble sort
void bubblesort(double a[], int array_size);

//quick sort for vector
int split(double a[], int low, int high);
void qsort(double a[], int low, int high );

/*summation for array*/
int sum(int* vector, int len);
double sum(double* vector, int len);

/*change the hist into point sets*/
double*  hist2points(double *hist, int diameter, int *numdata);

/*compute mean and variance for the given points*/
void meanvariance(int* dataPoints, int numdata, double* mean, double* var, int numlayer);

/*copy a subset from a give points*/
int* subcopy(double* dataPoints, int numdata, struct array1d reidx);
int* subcopy(double* dataPoints, int numdata, struct cellData reidx, int dim);
/*compute the distance between points*/
//double* funDist(double *pi1, double *pi2, int ninst, int col);
double* funDist(double *pinstances, int ninst, int col);
void funDist(double *pinstances, int ninst, int col, int* pdist);
/*compute the neighborhood for each point according to its IDs*/
struct array1d searchneighs(int ptid, double *pdist, int ninst, int radius, int pre_radius);
struct cellData searchneighs_ex(int r, int c, double *pinstances, int ninst,int dim, int radius, int numlayer);
//struct array1d searchneighs(int ptid, double *pdist, int ninst, int radius, int pre_radius, bool *visited);
struct array1d searchneighs(int ptid, int *pdist, int ninst, int radius, int pre_radius, int numlayer);
struct array1d searchneighs(int ptid, int *pdist, int ninst, int radius, int pre_radius, bool *visited);
/*search points in a given radius*/
struct array1d searchneighs(int r, int c, double *pinstances, int ninst,int dim, int radius, int numlayer);
struct array1d searchneighs(int r, int c, double *pinstances, int ninst,int dim, int radius, int flag, int numlayer);
struct array1d searchneighs2(int r, int c, double *pinstances, int ninst,int dim, int radius, int numlayer, double* radjust);

/*compute the difference between two sets*/
struct array1d setdiff(struct array1d idx, struct array1d idx2);

/*type to orientation*/
void type2orient(int type, int numOrient, int* ori1, int* ori2);

/*local maxima in the given region*/
double *local_maxima(double *h, int radius, int*retval, int*neighbors);

/*suppress local less than a given value*/ 
int* suppress(double *v, int length, double thresh);

/*neighborhood maximum*/
bool lmax(double **a, int i, int j, int rows, int cols, int nn);

/*write matrix to file*/
void writeArraytoFile(double *array, int row, int col, const char* filename);
void writeMatrixtoFile(double **array, int row, int col);

/*find = > < integers*/
int* find(int *olist, int olen, int oid, char s);

/*compute geometric mean*/
double geomean(double *vector, int n);

/*check overlap between lambda*/
bool eoverlap(struct cellData *ele, struct cellData *eleneigh, int dim);

#endif
