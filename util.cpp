#include "util.h"

/* Generating integer vector */
int *int_vector(int n)
{
    int *v; 
    v = (int*)malloc(n*sizeof(int));
    return v; 
}
/* Generating integer matrix */
int **int_matrix(int m, int n)
{
    int **mat; 
    int i; 
    mat = (int**)malloc(m*sizeof(int*)); 
    for (i=0; i<m; i++)
        mat[i] = int_vector(n); 
    return mat; 
}


double *double_vector(int n)
{
    double *v; 
    v = (double*)malloc(n*sizeof(double));
    return v; 
}
/* Generating integer matrix */
double **double_matrix(int m, int n)
{
    double **mat; 
    int i; 
    mat = (double**)malloc(m*sizeof(double*)); 
    for (i=0; i<m; i++)
        mat[i] = double_vector(n); 
    return mat; 
}



/* Free matrix space */
void free_matrix(void **mat, int m, int n)
{
        int i;
        for (i=0; i<m; i++)
              free(mat[i]);
        free(mat);
}

/*
// Best dilate by one solution
int** dilate(int** image){
    for (int i=0; i<image.length; i++){
        for (int j=0; j<image[i].length; j++){
            if (image[i][j] == 1){
                if (i>0 && image[i-1][j]==0) image[i-1][j] = 2;
                if (j>0 && image[i][j-1]==0) image[i][j-1] = 2;
                if (i+1<image.length && image[i+1][j]==0) image[i+1][j] = 2;
                if (j+1<image[i].length && image[i][j+1]==0) image[i][j+1] = 2;
            }
        }
    }
    for (int i=0; i<image.length; i++){
        for (int j=0; j<image[i].length; j++){
            if (image[i][j] == 2){
                image[i][j] = 1;
            }
        }
    }
    return image;
}
*/
// Best dilate by one solution
void imdilate(int *image, int imw, int imh){
    for (int i=0; i<imh; i++){
        for (int j=0; j<imw; j++){
            if (image[px2(i,j, imh,imw)] == 1){
                if (i>0 && image[px2(i-1, j, imh,imw)]==0) image[px2(i-1,j, imh,imw)] = 2;
                if (j>0 && image[px2(i,j-1,imh,imw)]==0) image[px2(i,j-1,imh,imw)] = 2;
                if (i+1<imh && image[px2(i+1,j,imh,imw)]==0) image[px2(i+1,j,imh,imw)] = 2;
                if (j+1<imw && image[px2(i,j+1,imh,imw)]==0) image[px2(i,j+1,imh,imw)] = 2;
            }
        }
    }
    for (int i=0; i<imw; i++){
        for (int j=0; j<imh; j++){
            if (image[px2(i,j,imw,imh)] == 2){
                image[px2(i,j,imw,imh)] = 1;
            }
        }
    }
    //return image;
}

double imoverlap(int *im1, int *im2, int imh, int imw){

int *imunion = (int*)malloc(imh*imw*sizeof(int));
for(int i=0;i<imh;i++){
for(int j=0; j< imw; j++)
{
    if(im1[px2(i,j,imw,imh)]>=im2[px2(i,j,imw,imh)])
    imunion[px2(i,j,imw,imh)] = im1[px2(i,j,imw,imh)];
    else
    imunion[px2(i,j,imw,imh)] =im2[px2(i,j,imw,imh)];
}
}

//compute the sum 
int sum1,sum2,sum;
sum1=0;
sum2=0;
sum= 0;
for(int i=0;i<imh;i++){
for(int j=0; j< imw; j++)
{
    sum1 = im1[px2(i,j,imw,imh)] + sum1;
    sum2 = im2[px2(i,j,imw,imh)] + sum2;
    sum = sum+ imunion[px2(i,j,imw,imh)];
}
}
return (sum1+sum2-sum)/(sum+0.001);
}



/*quick sort*/
void swap(double **a, int dim, int left, int right)
{
    double temp;
    for(int i=0; i< dim; i++)
    {   temp=a[left][i];
        a[left][i]=a[right][i];
        a[right][i]=temp;
    }
}//end swap
 
/*quick sort for matrix */
void quicksort(double **a,int dim, int low, int high)
{
     int pivot;
     // Termination condition! 
     if (high >low)
     {
        pivot = partition( a, dim, low, high );
        if(low<pivot-1)
        quicksort(a, dim, low, pivot-1);
        if(pivot+1<high)
        quicksort(a, dim, pivot+1, high);
     }
} //end quicksort
 
int partition( double **a,int dim, int low, int high )
{
     int left, right;
     double *pivot_item= (double*)malloc(dim*sizeof(double)); 
     int pivot = left = low; 
     
     for(int i=0; i<dim; i++)        
         pivot_item[i] = a[low][i]; 
     right = high;
     while ( left < right ) 
     {
      // Move left while item < pivot 
      while(a[left][0] <= pivot_item[0] && left<high) 
           left++;
      // Move right while item > pivot 
      while( a[right][0] > pivot_item[0] && right>low) 
           right--;
      if ( left < right ) 
            swap(a, dim, left,right);
     }
     // right is final position for the pivot    
     if(low<=high)
     for(int i=0; i< dim; i++){
        a[low][i] = a[right][i];
        a[right][i] = pivot_item[i];    
     }
     free(pivot_item);
     return right;
}//end partition


void bubblesort(double a[], int array_size)
{
     int i, j;
     double temp;
     for (i = 0; i < (array_size - 1); ++i)
     {
          for (j = 0; j < array_size - 1 - i; ++j )
          {
               if (a[j] > a[j+1])
               {
                    temp = a[j+1];
                    a[j+1] = a[j];
                    a[j] = temp;
               }
          }
     }
 }   


/*quick sort for array*/

int split(double a[], int low, int high)
{
   int pivot, left, right;
   double pivot_item, tmp;
   pivot_item = a[low];
   pivot = left = low;
   right = high;
   while ( left < right ){
     /* Move left while item < pivot */
     while( a[left] <= pivot_item && left<high) left++;
     /* Move right while item > pivot */
     while( a[right] > pivot_item && right>low) right--;
     if ( left < right ) 
         {tmp = a[left]; a[left] = a[right]; a[right] = tmp;}
   }
   /* right is final position for the pivot */
   if(right>0)
   {    
       a[low] = a[right];
       a[right] = pivot_item;
   }
   return right;
}

void qsort(double a[], int low, int high )
{
   int pivot;
   /* Termination condition! */
   if ( high > low )
   {
         pivot = split(a, low, high);
         if(low<pivot-1)
         qsort(a, low, pivot-1);
         if(pivot+1<high)
         qsort(a, pivot+1, high);
   }
}

/*good one works*/
void quickSort2(double arr[], int left, int right)
{
    int i = left, j = right;
    double tmp;
    double pivot = arr[abs((left + right) / 2)];
    /* partition */
    while (i < j) {
        while (arr[i] < pivot)
              i++;
        while (arr[j] > pivot)
              j--;
        if (i <= j) {
              tmp = arr[i];
              arr[i] = arr[j];
              arr[j] = tmp;
              i++;
              j--;
         }
    }

    /* recursion */
    if (left < j)
        quickSort2(arr, left, j);

    if (i< right)
        quickSort2(arr, i, right);
}



double sum(double* vector, int len)
{
   double ret = 0;
   for(int i=0; i<len; i++)
      ret = ret+ vector[i];
   return ret;
    
}


int sum(int* vector, int len)
{
   int ret = 0;
   for(int i=0; i< len; i++)
      ret = ret+ vector[i];
   return ret;
    
}


/*change the hist into point sets*/
double*  hist2points(double *hist, int diameter, int *numdata)
{
    *numdata = sum(hist, diameter*diameter);
    int dim = 2;
    double *dataPoints= (double*)malloc((*numdata)*dim*sizeof(double));
    int idx = 0;
    for(int i=0; i < diameter; i++)
    for(int j=0; j < diameter; j++)
    {
        for(int k=1; k<=hist[px2(i, j, diameter, diameter)]; k++)
        //for(int k=1; k<=hist[i*diameter + j]; k++)   
        {    
            //dataPoints[idx*2+0] = i+1;
            //dataPoints[idx*2+1] = j+1;
            
            dataPoints[px2(idx,0, *numdata,dim)] = i+1;
            dataPoints[px2(idx,1, *numdata,dim)] = j+1;
            
            
            idx = idx + 1;
        }
    }
    return dataPoints;
}


/*copy a subset from a give points*/
int* subcopy(double* dataPoints, int numdata, struct array1d reidx)
{
    int dim = 2;
    int* subdataPoints= (int*)malloc(dim*reidx.length*sizeof(int));
    for(int i=0; i < reidx.length; i++)
    {
     
        subdataPoints[px2(i,0,reidx.length,dim)] = (int)dataPoints[px2(reidx.a[i],0,numdata,dim)];
        subdataPoints[px2(i,1,reidx.length,dim)] = (int)dataPoints[px2(reidx.a[i],1,numdata,dim)];
    }
    
    return subdataPoints;
}

/*copy a subset from a give points*/
int* subcopy(double* dataPoints, int numdata, struct cellData reidx, int dim)
{
    if(dim != 2)
        printf("error! the dimension should be 2!\n");
    int* subdataPoints= (int*)malloc(dim*reidx.length*sizeof(int));
    for(int i=0; i < reidx.length; i++)
    {
     
        subdataPoints[px2(i,0,reidx.length,dim)] = (int)dataPoints[px2((int)reidx.a[i],0,numdata,dim)];
        subdataPoints[px2(i,1,reidx.length,dim)] = (int)dataPoints[px2((int)reidx.a[i],1,numdata,dim)];
    }
    
    return subdataPoints;
}

/*compute mean and variance for the given points*/
void meanvariance(int* dataPoints, int numdata, double* mean, double* var, int numlayer)
{
    int dim = 2;
    /*mean*/
    mean[0] = 0; mean[1] =0;
    var[0] = 0; var[1] = 0;
    for(int i=0; i< numdata; i++)
    {
        mean[0] = mean[0] + dataPoints[px2(i,0,numdata,dim)];
        mean[1] = mean[1] + dataPoints[px2(i,1,numdata,dim)];
    
    }
    mean[0]= mean[0]/numdata;mean[1]= mean[1]/numdata;
    /*variance here*/
    for(int i=0; i< numdata; i++)
    {
        var[0] = var[0] + (dataPoints[px2(i,0,numdata,dim)]- mean[0])*(dataPoints[px2(i,0,numdata,dim)]- mean[0]);
        var[1] = var[1] + (dataPoints[px2(i,1,numdata,dim)]- mean[1])*(dataPoints[px2(i,1,numdata,dim)]- mean[1]);
    
    }
    if(numdata>2)
    {
        var[0]= var[0]/(numdata-1)+1.5;
        var[1]= var[1]/(numdata-1)+1.5;
    }
    else{
        var[0]= 1.5*numlayer;var[1]= 1.5*numlayer;
    }
    
}

/*compute the distance between points*/
double* funDist(double *pinstances, int ninst, int col)
{
    double* pdist = (double*)malloc(ninst*ninst*sizeof(double));
    int i,j;
    int dx, dy;
    #pragma omp parallel for 
    for(i=0; i< ninst; i++){
        for(j=i+1;j<ninst; j++)
        {
            dy = (pinstances[px2(i, 1, ninst, col)]-pinstances[px2(j, 1, ninst, col)]);
            dx = (pinstances[px2(i, 2, ninst, col)]-pinstances[px2(j, 2, ninst, col)]);
            pdist[px2(i,j, ninst, ninst)] = sqrt(dx*dx + dy*dy);
            pdist[px2(j,i, ninst, ninst)] = pdist[px2(i,j, ninst, ninst)];
        }
        pdist[px2(i,i, ninst, ninst)] = 0;
    }
    //pdist[px2(i,i, ninst, ninst)] = 0;
    return pdist;
}



void funDist(double *pinstances, int ninst, int col, int* pdist)
{

    int i,j;
    int dx, dy;
    #pragma omp parallel for 
    for(i=0; i< ninst; i++){
        for(j=i+1;j<ninst; j++)
        {
            dy = (pinstances[px2(i, 1, ninst, col)]-pinstances[px2(j, 1, ninst, col)]);
            dx = (pinstances[px2(i, 2, ninst, col)]-pinstances[px2(j, 2, ninst, col)]);
            pdist[px2(i,j, ninst, ninst)] = round(sqrt(dx*dx + dy*dy));
            pdist[px2(j,i, ninst, ninst)] = pdist[px2(i,j, ninst, ninst)];
        }
        pdist[px2(i,i, ninst, ninst)] = 0;
    }
    
}



struct array1d searchneighs(int ptid, double *pdist, int ninst, int radius, int pre_radius)
{   
    int num=0;
    struct array1d idx;
    idx.length = -1;
    for(int iinst = 0; iinst<ninst; iinst++)
    {
          if((pdist[px2(ptid, iinst, ninst, ninst)] <radius) && (pdist[px2(ptid, iinst, ninst, ninst)] >=pre_radius))
          {

            if(idx.length >MAX_LEN-2)
                break;
             idx.length= idx.length+1;
             idx.a[idx.length] = iinst;

          }

    }
 
   idx.length= idx.length+1;
   return idx;
}

struct array1d searchneighs(int ptid, int *pdist, int ninst, int radius, int pre_radius, int numlayer)
{   
    int num=0;
    struct array1d idx;
    idx.length = -1;
    for(int iinst = 0; iinst<ninst; iinst++)
    {
          if((pdist[px2(ptid, iinst, ninst, ninst)] <radius) && (pdist[px2(ptid, iinst, ninst, ninst)] >=pre_radius))
          {

            if(idx.length >MAX_LEN-2)
                break;
             idx.length= idx.length+1;
             idx.a[idx.length] = iinst;

          }

    }
 
   idx.length= idx.length+1;
   return idx;
}


struct array1d searchneighs(int ptid, int *pdist, int ninst, int radius, int pre_radius, bool *visited)
{   
    int num=0;
    struct array1d idx;
    idx.length = -1;
    for(int iinst = 0; iinst<ninst; iinst++)
    {
        if(!visited[iinst] && ptid!=iinst){
          if((pdist[px2(ptid, iinst, ninst, ninst)] <radius) && (pdist[px2(ptid, iinst, ninst, ninst)] >=pre_radius))
          {

            if(idx.length >MAX_LEN-2)
                break;
             idx.length= idx.length+1;
             idx.a[idx.length] = iinst;

          }
        }
    }
 
   idx.length= idx.length+1;
   return idx;
}




/*search points in a given radius & store index*/
struct array1d searchneighs(int r, int c, double *pinstances, int ninst,int dim, int radius, int numlayer)
{   
    int num=0;
    struct array1d idx;
    idx.length = -1;
    if(dim >=3) 
    {
        for(int iinst = 0; iinst<ninst; iinst++)
        {
             /*if((pinstances[px2(iinst, 1, ninst, dim)]>= (r-radius+1)) && (pinstances[px2(iinst, 1, ninst, dim)]<= (r+radius-1))
                     && pinstances[px2(iinst, 2, ninst, dim)]>= (c-radius+1) && pinstances[px2(iinst, 2, ninst, dim)]<= (c+radius-1))*/
              if((abs(pinstances[px2(iinst, 1, ninst, dim)]-r) <radius) && (abs(pinstances[px2(iinst, 2, ninst, dim)]-c) <radius))
              {

                if(idx.length >MAX_LEN-2)
                    break;
                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;

              }

         }
    }
    else
    {
        
        for(int iinst = 0; iinst<ninst; iinst++)
        {
             if((pinstances[px2(iinst, 0, ninst, dim)]>= (r-radius+1)) && (pinstances[px2(iinst, 0, ninst, dim)]<= (r+radius-1))
                     && (pinstances[px2(iinst, 1, ninst, dim)]>= (c-radius+1)) && (pinstances[px2(iinst, 1, ninst, dim)]<= (c+radius-1)))
              {

                if(idx.length >MAX_LEN-2)
                    break;
                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;

              }

         }
        
        
    }
   idx.length= idx.length+1;
   return idx;
}



/*search points in a given radius & store index*/
struct array1d searchneighs(int r, int c, double *pinstances, int ninst,int dim, int radius, int flag, int numlayer)
{   
    int num=0;
    struct array1d idx;
    idx.length = -1;
    if(dim >=3 && flag) 
    {
        for(int iinst = 0; iinst<ninst; iinst++)
        {
             /*if((pinstances[px2(iinst, 1, ninst, dim)]>= (r-radius+1)) && (pinstances[px2(iinst, 1, ninst, dim)]<= (r+radius-1))
                     && pinstances[px2(iinst, 2, ninst, dim)]>= (c-radius+1) && pinstances[px2(iinst, 2, ninst, dim)]<= (c+radius-1))*/
              if((abs(pinstances[px2(iinst, 1, ninst, dim)]-r) <radius) && (abs(pinstances[px2(iinst, 2, ninst, dim)]-c) <radius) && pinstances[px2(iinst, 2, ninst, dim)]==1)
              {

                if(idx.length >MAX_LEN-2)
                    break;
                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;

              }

         }
    }
    else if(dim >=3) 
    {
        for(int iinst = 0; iinst<ninst; iinst++)
        {

              if((abs(pinstances[px2(iinst, 1, ninst, dim)]-r) <radius) && (abs(pinstances[px2(iinst, 2, ninst, dim)]-c) <radius))
              {

                if(idx.length >MAX_LEN-2)
                    break;
                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;

              }

         }
    }
    else
    {
        
        for(int iinst = 0; iinst<ninst; iinst++)
        {
             if((pinstances[px2(iinst, 0, ninst, dim)]>= (r-radius+1)) && (pinstances[px2(iinst, 0, ninst, dim)]<= (r+radius-1))
                     && (pinstances[px2(iinst, 1, ninst, dim)]>= (c-radius+1)) && (pinstances[px2(iinst, 1, ninst, dim)]<= (c+radius-1)))
              {

                if(idx.length >MAX_LEN-2)
                    break;
                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;

              }

         }
        
        
    }
   idx.length= idx.length+1;
   return idx;
}


/*handle any number of points search points in a given radius & store index*/
struct cellData searchneighs_ex(int r, int c, double *pinstances, int ninst,int dim, int radius, int numlayer)
{
    int num=0;
    struct cellData idx;
    idx.length = -1;
    idx.a = (double*)malloc(ninst*sizeof(double));
    if(dim >=3)
    {
        for(int iinst = 0; iinst<ninst; iinst++)
        {
             /*if((pinstances[px2(iinst, 1, ninst, dim)]>= (r-radius+1)) && (pinstances[px2(iinst, 1, ninst, dim)]<= (r+radius-1))
                     && pinstances[px2(iinst, 2, ninst, dim)]>= (c-radius+1) && pinstances[px2(iinst, 2, ninst, dim)]<= (c+radius-1))*/
              if((abs(pinstances[px2(iinst, 1, ninst, dim)]-r) <radius) && (abs(pinstances[px2(iinst, 2, ninst, dim)]-c) <radius))
              {

                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;

              }

         }
    }
    else
    {

        for(int iinst = 0; iinst<ninst; iinst++)
        {
             if((pinstances[px2(iinst, 0, ninst, dim)]>= (r-radius+1)) && (pinstances[px2(iinst, 0, ninst, dim)]<= (r+radius-1))
                     && (pinstances[px2(iinst, 1, ninst, dim)]>= (c-radius+1)) && (pinstances[px2(iinst, 1, ninst, dim)]<= (c+radius-1)))
              {

                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;

              }

         }


    }
   idx.length= idx.length+1; 
   return idx;
}                                   


/*search points in a given radius & store index*/
struct array1d searchneighs2(int r, int c, double *pinstances, int ninst,int dim, int radius, int numlayer, double* radjust)
{   
    int num=0;
    struct array1d idx;
    idx.length = -1;
    //adust the radius here
    radius = radius - radjust[numlayer];
    //idx = malloc(MAX_LEN*sizeof(int));
    if(dim >=3) //more than three dimenstion structure
    {
        for(int iinst = 0; iinst<ninst; iinst++)
        {
             //if (abs(r- (pinstances[px2(iinst, 1, ninst, dim)]-radjust[numlayer]) + 1) <= radius && abs(c- (pinstances[px2(iinst, 2, ninst, dim)]-radjust[numlayer]) + 1) <= radius)
            if((pinstances[px2(iinst, 1, ninst, dim)]>= r-radius+1) && (pinstances[px2(iinst, 1, ninst, dim)]<= r+radius-1)
                     && (pinstances[px2(iinst, 2, ninst, dim)]>= c-radius+1) && (pinstances[px2(iinst, 2, ninst, dim)]<= c+radius-1))
            {
                 
                if(idx.length >MAX_LEN-1)
                    break;
                 
                 
                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;
                
             }

         }
    }
    else // two dimention points
    {  
        for(int iinst = 0; iinst<ninst; iinst++)
        {
             //if (abs(r- (pinstances[px2(iinst, 0, ninst,dim)]-radjust[numlayer])  +1) <= radius && abs(c- (pinstances[px2(iinst, 1, ninst,dim)]-radjust[numlayer])  +1) <= radius)
             if((pinstances[px2(iinst, 0, ninst, dim)]>= r-radius+1) && (pinstances[px2(iinst, 0, ninst, dim)]<= r+radius-1)
                     && (pinstances[px2(iinst, 1, ninst, dim)]>= c-radius+1) && (pinstances[px2(iinst, 1, ninst, dim)]<= c+radius-1))
            {
                 if(idx.length >MAX_LEN-1)
                    break;
                 
                 idx.length= idx.length+1;
                 idx.a[idx.length] = iinst;     
                   
             }

         } 
        
    }
    idx.length= idx.length+1;
    return idx;
}


/*get the difference between two sets*/
struct array1d setdiff(struct array1d idx, struct array1d idx2)
{
    int *indicator = (int*)malloc(sizeof(int)*idx.length);
    int i,j;
    for(i=0; i<idx.length; i++)
    {
        indicator[i] = 1;
    }
    
    for(i=0; i<idx.length; i++)
    {

        for(j=0;j<idx2.length; j++)
        {
              if(idx.a[i]== idx2.a[j])
                    indicator[i] = 0;
        }

     }

     int numel = -1;
     struct array1d idxdiff;
     idxdiff.length = 0;
     for(i=0; i<idx.length; i++)
     {
         if(indicator[i]==1)
         {
            numel = numel + indicator[i];
            idxdiff.a[numel] = idx.a[i];
         }
     }
     /*
     while(i<idx.length)
     {
        numel = numel + indicator[i];
        idxdiff.a[numel] = idx.a[i];
        i = i+1;
     }
     */
     idxdiff.length = numel+1;
     free(indicator);
     return idxdiff;
}

/*any overlap between two sets*/
bool eoverlap(struct cellData *ele, struct cellData *eleneigh, int dim){
    
    int i, j;
    bool found = false;
    for(i=0; i< ele->length; i++)
        for(j=0; j<eleneigh->length; j++)
        {
            if(ele->a[i*dim+0] == eleneigh->a[j*dim+0] && ele->a[i*dim+1] == eleneigh->a[j*dim+1])
            {
                found = true;
            	return found;
            }
        }
    return found;
}


void type2orient(int type, int numOrient, int* ori1, int* ori2)
{
    //ori2 = mod(type, numOrient);
    /*
    *ori2 = (type-1)%numOrient +1;
    *ori1 = type - (*ori1-1)*numOrient;
    */
    *ori1 = floor((type-1)/numOrient) +1;
    *ori2 = type - (*ori1-1)*numOrient;
    
    
}




double *local_maxima(double* h, int radius, int* retval, int* neighbors)
{
     //int MAX_POINTS = 100;
     int i,j;
     double sigma = 2*radius*0.07;
     double **gKernel;
     int klen = round((5*sigma-1)/2);//radius of the kernel
     
     createFilter(klen, sigma, gKernel);
     /*
     double **a;
     conv(h, 2*radius-1, 2*radius-1, gKernel,2*klen+1, 2*klen+1, a);
      */
     double *a = conv2(h, 2*radius-1, 2*radius-1, gKernel,2*klen+1, 2*klen+1);
     //writeArraytoFile(a, 2*radius-1, 2*radius-1, "conv2.txt");
     // find the maximum in the histogram
     double maximum = 0;
     for(int i=0;i<2*radius-1;i++)
         for(int j=0; j<2*radius-1; j++)
             /*if(a[i][j] >maximum)
                 maximum = a[i][j];*/
             if(a[px2(i,j, 2*radius-1, 2*radius-1)]>maximum)
                 maximum = a[px2(i,j, 2*radius-1, 2*radius-1)];
     // find the potentional points
     int dim = 3;
     double* points = (double*)malloc(MAX_POINTS*dim*sizeof(double));
     
     int count = 0;
     for(i=0;i<2*radius-1;i++)
         for(j=0; j<2*radius-1; j++)
         {
             //if(a[i][j]>0.05*maximum)
             if(a[px2(i,j, 2*radius-1, 2*radius-1)]>0.05*maximum)
             {   /* 
                 if((i-1)>=0 && j-1>=0 && a[i-1][j-1]<= a[i][j] && (i-1)>=0 && a[i-1][j]<= a[i][j] 
                         && (i-1)>=0 && (j+1)< 2*radius-1 && a[i-1][j+1]<= a[i][j] && a[i][j-1]<= a[i][j] && a[i][j+1] <= a[i][j] 
                         && (i+1)<2*radius && a[i+1][j] <= a[i][j] && a[i+1][j-1] <= a[i][j] && a[i+1][j+1]<=a[i][j])
                  */
                 // processing neighbors
                 /*if(lmax(a, i, j, 2*radius-1, 2*radius-1, 4))*/
                int in, jn;
                int nnc = 0;
                int nn =4;
                for(int k=0; k< nn;k++)
                {

                    in = MIN(MAX(neighbors[k*2+0]+i,0), 2*radius-2);
                    jn = MIN(MAX(neighbors[k*2+1]+j,0), 2*radius-2);
                    //if(a[i][j]> a[in][jn])
                    if(a[px2(i,j, 2*radius-1, 2*radius-1)]>a[px2(in,jn, 2*radius-1, 2*radius-1)])
                        nnc++;

                }
                 
                 if(nnc==nn)
                 {
                     if (count > MAX_POINTS-1)
                         break;
                     else
                     {
                         points[count*dim+0] = i+1;
                         points[count*dim+1] = j+1;
                         points[count*dim+2] = a[px2(i,j, 2*radius-1, 2*radius-1)]; //a[i][j];
                         count++;
                     }
                 }
             }
         }
     //free allocated matrix
     free_matrix((void**)gKernel, 2*klen+1, 2*klen+1);
     // free_matrix((void**)a, 2*radius-1, 2*radius-1);
     free(a);
     //sort the result
      double**final = (double**)malloc(count*sizeof(double*));
      for(i=0; i<count;i++)
          final[i] = (double*)malloc(3*sizeof(double));
      for(i=0; i< count; i++)
      {   
          final[i][0] = -points[i*dim+2];
          final[i][1] = points[i*dim+0];
          final[i][2] = points[i*dim+1];
      }
      /*free allocation*/
      if(points) {free(points);points = NULL;}
      /*sort according to the score*/
      quicksort(final,3, 0, count-1);
      /* non maximum suppression*/
      double thresh = (2*radius-1)*0.24;
      /*compute the distance between two points*/
      double *dist = (double*)malloc(count*count*sizeof(double));
      for(i=0; i< count; i++)
      {    
          for(int j=0; j < count; j++)
          {
               dist[i*count+j] = sqrt((final[i][1]-final[j][1])*(final[i][1]-final[j][1]) + (final[i][2]-final[j][2])*(final[i][2]-final[j][2]));
               
          
          }
          dist[i*count+i] = INF;    
      }
      int *indicator = (int*)malloc(count*sizeof(int));
      for(i=0; i<count; i++)
          indicator[i] = 1;
      /**delete too close points**/
      i = 0;
      do
      {
          
          for(j=i+1; j<count; j++)
          {
              if (dist[i*count +j] < thresh)
              {
                  indicator[j] = 0;
                  dist[i*count+j] = INF;
              }

          }
          i = i+1;
          
      }while(i<count);
      /*return the final results*/
      *retval = sum(indicator, count);
      double*peak = (double*)malloc((*retval)*dim*sizeof(double)); 
      // copy it to peak
      int idx = 0;
      for(i= 0; i< count; i++)
      {    
          if(indicator[i] ==1)
          {
            peak[idx*dim+0] = -final[i][0];
            peak[idx*dim+1] = final[i][1];
            peak[idx*dim+2] = final[i][2];
            idx = idx+1;
          }
      }
      // release memory
      
      free(dist);
      free(final);
      free(indicator);
      return peak;
}
 
 

int* suppress(double *v, int length, double thresh)
{
    int *indicator = (int*)malloc(sizeof(int)*length);
    int num = 0;
    for(int i=0; i<length; i++)
    {
        if(v[i] < thresh)
            indicator[num++] = i;
    
    }
    return indicator;
}

bool lmax(double **a, int i, int j, int rows, int cols, int nn)
{
    int **neighbors = int_matrix(nn,2);
    if(nn==4) //4-neighborhood
    {   //neighbors = {{-1,0}, {0, -1}, {0,+1}, {+1,0}};
    
        neighbors[0][0] = -1; neighbors[0][1] = 0; 
        neighbors[1][0] = 0; neighbors[1][1] = -1; 
        neighbors[2][0] = 0; neighbors[2][1] = 1; 
        neighbors[3][0] = 1; neighbors[3][1] = 0; 
    }
    else if(nn==8)
    {//8 - neighborhood
     // neighbors = {{-1,-1}, {-1,0}, {-1,+1}, {0, -1}, {0,+1}, {+1,-1},{+1,0},{+1,+1}};
        
        neighbors[0][0] = -1; neighbors[0][1] = 0; 
        neighbors[1][0] = 0; neighbors[1][1] = -1; 
        neighbors[2][0] = 0; neighbors[2][1] = 1; 
        neighbors[3][0] = 1; neighbors[3][1] = 0; 
        neighbors[4][0] = -1; neighbors[4][1] = -1; 
        neighbors[5][0] = -1; neighbors[5][1] = 1; 
        neighbors[6][0] = 1; neighbors[6][1] = 1; 
        neighbors[7][0] = 1; neighbors[7][1] = -1; 
        
        
    }
    else
    {
        printf("specify only 4 0r 8 neighs");
        exit(0);
    }
    int in;
    int jn;
    int count = 0;
    for(int k=0; k< nn;k++)
    {
    
        in = MIN(MAX(neighbors[k][0]+i,0), rows-1);
        jn = MIN(MAX(neighbors[k][1]+j,0), cols-1);
        if(a[i][j]>=a[in][jn])
            count++;
    
    }
    free_matrix((void**)neighbors, nn,2);
    if(count==nn)
        return true;
    else
        return false;
}



void writeArraytoFile(double *array, int row, int col, const char *filename)
{
    if(filename==NULL)
        return;
    FILE *matrix=fopen(filename, "w");
    if(matrix==NULL)
        return;
    int i=0, j=0;
    for(i=0;i<row;i++) 
    {
        for(j=0;j<col;j++) 
        { 
            fprintf(matrix, "%lf\t", array[px2(i,j,row,col)]);
        }
        fprintf(matrix, "\n");
    }

    fclose(matrix);
    
    /*
    for(j= 0; j<row; j++)
    {
        for(k=0; k<col; k++)
        {    
            printf("%d ", array[px2(j, k, diameter, diameter)]);
        }
        printf("\n");
    }
    */
}


void writeMatrixtoFile(double **array, int row, int col)
{
    FILE *matrix=fopen("matrix_cov.txt", "w");
    int i=0, j=0;
    for(i=0;i<row;i++) 
    {
        for(j=0;j<col;j++) 
        { 
            fprintf(matrix, "%lf\t", array[i][j]);
        }
        fprintf(matrix, "\n");
    }

    fclose(matrix);
}


/*find elements that equal, larger or small then a given id*/
int* find(int *olist, int olen, int oid, char s)
{
    int *indicator = (int*)malloc(olen*sizeof(int));
    int i;
    for(i = 0; i< olen; i++)
    { 
        indicator[i] = 0;
    }
    if (s=='e') //find equal
    {
         for(i= 0; i< olen; i++)
         {
             if(olist[i] == oid)
                 indicator[i] = 1;
         }
    }
    else if(s=='l') //larger
    {
         for(i= 0; i< olen; i++)
         {
             if(olist[i] > oid)
                 indicator[i] = 1;
         }
    }
    else if(s=='s') //small
    {
         for(i= 0; i< olen; i++)
         {
             if(olist[i] < oid)
                 indicator[i] = 1;
         }
    }
    return indicator;
}


double geomean(double *vector, int n)
{   
    double value =1;
    for(int i=0; i<n; i++)
        value = value*vector[i];
    
    value = pow(value, 1/(double)n);
    return value;
    
}




/*
int* nonmaxsuppress(double *v, int length, double thresh)
{
    int *indicator = (int*)malloc(sizeof(int)*length);
    int num = 0;
    for(int i=0; i<length; i++)
    {
        if(v[i] < thresh)
            indicator[num++] = i;
    
    }
    return indicator;
}
*/
