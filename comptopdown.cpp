#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "comptopdown.h"
/*define the total number of layers*/
// int totLayers = 7;
// int numsub = 2;
//#define px2(x, y, lengthx, lengthy) ((x) + (y)*(lengthx) - 0) 
//#define px3(x, y, z, lengthx, lengthy, lengthz) (z+ ((x) + ((y))*(lengthx))*lengthz - 0) 
// struct Parts model[totLayers];
/*single image inference*/
Data* comptopdown(double *morients, int imh, int imw, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust)
{
    
     Data *output = (Data*)malloc(sizeof(Data));
     int i,j,k;
     int id;
     int ninst, col, iinst, newid, ncount;
     int r, c, o1, o2;
     struct array1d idx, idx1;
     int pDim = 4;
     ninst = 0;
     /*get the number of sampled points*/
     for(i=0; i< imh; i++)
         for(j=0;j<imw;j++)
             if(morients[px2(i,j,imh,imw)] >0)
             ninst++;
     /*construct pinstances*/ 
     double *pinstances = (double*)malloc(pDim*ninst*(sizeof(double)));
     id = 0; 
     for(i=0; i< imh; i++)
         for(j=0;j<imw;j++){
             if(morients[px2(i,j,imh,imw)] >0)
             {
             pinstances[px2(id, 0, ninst, pDim)] = morients[px2(i,j, imh, imw)];
             pinstances[px2(id, 1, ninst, pDim)] = i+1;//y coordinate
             pinstances[px2(id, 2, ninst, pDim)] = j+1;//x coordinate
             pinstances[px2(id, 3, ninst, pDim)] = 1;//score
             id = id+1;}
     }

    /*define the output and allocate space*/
    int dimsOutput = 4;
    // allocate memory
    output[0].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
    output[0].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

    if(ninst==0)
    {
        output[0].pinstances[0].row = 0;
        output[0].pinstances[0].col = 0;
        output[0].lambda[0].numel =0;
        return output;
    }
    
    int radius = model[numlayer-1].radius;//get the radius from the model
    int pre_radius = 0;
    bool *indicator = (bool*)malloc(ninst*sizeof(bool));
    for(iinst =0; iinst < ninst; iinst++)
        indicator[iinst] = false;
    struct Node *point = (struct Node*)malloc(ninst*sizeof(struct Node));
    /*counting the matched shaped*/
    ncount = 0;
    for(iinst = 0; iinst< ninst; iinst++) //match model by window scanning
    {
        r = pinstances[px2(iinst, 1, ninst, pDim)];
        c = pinstances[px2(iinst, 2, ninst, pDim)];
        o1 = pinstances[px2(iinst, 0, ninst, pDim)];
         // find the neighbors 
        idx = searchneighs(r, c, pinstances, ninst, pDim, radius, numlayer);
        pre_radius = 2;
        if(numlayer>1)
            pre_radius = model[numlayer-2].radius;
        idx1 = searchneighs2(r, c, pinstances, ninst, pDim, pre_radius, numlayer-1, radjust);
        // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
        idx = setdiff(idx,idx1);
        
        if(idx.length>0)
        {
            int *y = (int*)malloc(sizeof(int)*idx.length); 
            //int x[idx.length];
            int *x = (int*)malloc(sizeof(int)*idx.length); 
            int *iy = (int*)malloc(sizeof(int)*idx.length); 
            int *ix = (int*)malloc(sizeof(int)*idx.length); 
            int *o2list = (int*)malloc(sizeof(int)*idx.length); 
            int *dyx = (int*)malloc(2*idx.length*sizeof(int));
            for(i=0;i< idx.length; i++)
            {

                y[i] = (int)pinstances[px2(idx.a[i], 1, ninst, pDim)];//pinstances[idx.a[i]][2];
                x[i] = (int)pinstances[px2(idx.a[i], 2, ninst, pDim)];//pinstances[idx.a[i]][3]; 
                o2list[i] = (int)pinstances[px2(idx.a[i], 0, ninst, pDim)];//pinstances[idx.a[i]][1];
                iy[i] = y[i] - r + radius;
                ix[i] = x[i] -c + radius;
                dyx[2*i]= iy[i]- radius;
                dyx[2*i+1] = ix[i]-radius; 

            }

            int numdyx = idx.length;
            printf("check the point id: %d\n", iinst);
            bool found = topdownmatch(pinstances, ninst, pDim, o1, r, c, dyx, o2list, numdyx, model, numlayer, numOrient, 
                    numsub, radjust, discount, &point[iinst]);

            indicator[iinst] = found;
            if(found)     /*counting in the inference stage*/
                ncount++;
            if(y!=NULL)
            {   free(y); y=NULL;}
            if(x!=NULL)
            {   free(x); x= NULL;}
            if(iy!=NULL)
            {   free(iy); iy= NULL;}
            if(ix!= NULL)
            {   free(ix); ix = NULL;}
            if(o2list!= NULL)
            {   free(o2list); o2list = NULL;}
            if(dyx!= NULL)
            {   free(dyx);dyx = NULL;}

        }      

    }
    
    /*new pid in the inference stage*/
    newid = 0;
    // allocate memory for elements
    output[0].pinstances[0].elements = (double*)malloc(ncount*dimsOutput*sizeof(double));
    output[0].lambda[0].elements  =  (struct cellData*)malloc(ncount*sizeof(struct cellData));
    col = dimsOutput;
    // copy results from the point
    for(iinst =0; iinst < ninst; iinst++)
        if(indicator[iinst]){ //update a match
                
                output[0].pinstances[0].elements[px2(newid, 3, ncount, col)] = point[iinst].score; //matched score
                output[0].pinstances[0].elements[px2(newid, 0, ncount, col)] = point[iinst].compid;// matched top level partid
                output[0].pinstances[0].elements[px2(newid, 1, ncount, col)] = pinstances[px2(iinst, 1, ninst, col)];//y
                output[0].pinstances[0].elements[px2(newid, 2, ncount, col)] = pinstances[px2(iinst, 2, ninst, col)];//x
                //allocate memory for its children
                output[0].lambda[0].elements[newid].a = (double*)malloc(2*point[iinst].length*sizeof(double));// y and x: two dimensions
                for(j=0; j< 2*point[iinst].length; j++)
                    output[0].lambda[0].elements[newid].a[j] =  point[iinst].a[j];
                output[0].lambda[0].elements[newid].length = point[iinst].length;
                newid = newid+1;
                
    }
    
    output[0].pinstances[0].row = ncount;
    output[0].pinstances[0].col = dimsOutput;
    output[0].lambda[0].numel = ncount;
    if(indicator!=NULL)
    {  
        free(indicator); 
        indicator = NULL;
    }
    if(point!=NULL)
    {    free(point); point = NULL; }
    
    if(pinstances!=NULL)
    {   free(pinstances);
        pinstances = NULL;
    }
    
    return output;
}

/*inference on same image with different scales with considering rotation for match*/
void comptopdown(Data *output, unsigned char* morients, int imh, int imw, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, double rotation)
{
    
     //Data *output = (Data*)malloc(sizeof(Data));
     int i,j,k;
     int id;
     int ninst, col, iinst, newid, ncount;
     int r, c, o1, o2;
     struct array1d idx, idx1;
     int pDim = 4;
     ninst = 0;
     /*get the number of sampled points*/
     for(i=0; i< imh; i++)
         for(j=0;j<imw;j++)
             if(morients[px2(i,j,imh,imw)] >0)
             ninst++;
     /*construct pinstances*/ 
     double *pinstances = (double*)malloc(pDim*ninst*(sizeof(double)));
     id = 0; 
     for(i=0; i< imh; i++)
         for(j=0;j<imw;j++){
             if(morients[px2(i,j,imh,imw)] >0)
             {
             pinstances[px2(id, 0, ninst, pDim)] = morients[px2(i,j, imh, imw)];
             pinstances[px2(id, 1, ninst, pDim)] = i+1;//y coordinate
             pinstances[px2(id, 2, ninst, pDim)] = j+1;//x coordinate
             pinstances[px2(id, 3, ninst, pDim)] = 1;//score
             id = id+1;}
     }

    /*define the output and allocate space*/
    int dimsOutput = 4;
    // allocate memory
    output[0].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
    output[0].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

    if(ninst==0)
    {
        output[0].pinstances[0].row = 0;
        output[0].pinstances[0].col = 0;
        output[0].lambda[0].numel =0;
        return;
    }
    
    int radius = model[numlayer-1].radius;//get the radius from the model
    int pre_radius = 0;
    bool *indicator = (bool*)malloc(ninst*sizeof(bool));
    for(iinst =0; iinst < ninst; iinst++)
        indicator[iinst] = false;
    struct Node *point = (struct Node*)malloc(ninst*sizeof(struct Node));
    /*counting the matched shaped*/
    ncount = 0;
    double *pdist = funDist(pinstances, ninst, pDim);
    #pragma omp parallel for 
    for(iinst = 0; iinst< ninst; iinst++) //match model by window scanning
    {
        r = pinstances[px2(iinst, 1, ninst, pDim)];
        c = pinstances[px2(iinst, 2, ninst, pDim)];
        o1 = pinstances[px2(iinst, 0, ninst, pDim)];
        //delete selected points by setting flag
        pinstances[px2(iinst, 3, ninst, pDim)] =0;
         // find the neighbors 
        /*idx = searchneighs(r, c, pinstances, ninst, pDim, radius, numlayer);
        pre_radius = 2;
        if(numlayer>1)
            pre_radius = model[numlayer-2].radius;
        idx1 = searchneighs2(r, c, pinstances, ninst, pDim, pre_radius, numlayer-1, radjust);
        // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
        idx = setdiff(idx,idx1);
        */
        pre_radius = 4;
        if(numlayer>1)
            pre_radius = model[numlayer-2].radius - radjust[numlayer-2];
        idx = searchneighs(iinst, pdist, ninst, radius, pre_radius);
        if(idx.length>0)
        {
            int *y = (int*)malloc(sizeof(int)*idx.length); 
            //int x[idx.length];
            int *x = (int*)malloc(sizeof(int)*idx.length); 
            int *iy = (int*)malloc(sizeof(int)*idx.length); 
            int *ix = (int*)malloc(sizeof(int)*idx.length); 
            int *o2list = (int*)malloc(sizeof(int)*idx.length); 
            int *dyx = (int*)malloc(2*idx.length*sizeof(int));
            for(i=0;i< idx.length; i++)
            {

                y[i] = (int)pinstances[px2(idx.a[i], 1, ninst, pDim)];//pinstances[idx.a[i]][2];
                x[i] = (int)pinstances[px2(idx.a[i], 2, ninst, pDim)];//pinstances[idx.a[i]][3]; 
                o2list[i] = (int)pinstances[px2(idx.a[i], 0, ninst, pDim)];//pinstances[idx.a[i]][1];
                iy[i] = y[i] - r + radius;
                ix[i] = x[i] -c + radius;
                dyx[2*i]= iy[i]- radius;
                dyx[2*i+1] = ix[i]-radius; 

            }

            int numdyx = idx.length;
            
            bool found = topdownmatch(pinstances, ninst, pDim, o1, r, c, dyx, o2list, numdyx, model, numlayer, numOrient, 
                    numsub, radjust, discount,rotation, &point[iinst]);


            indicator[iinst] = found;
            if(found)     /*counting in the inference stage*/
                ncount++;
            if(y!=NULL)
            {   free(y); y=NULL;}
            if(x!=NULL)
            {   free(x); x= NULL;}
            if(iy!=NULL)
            {   free(iy); iy= NULL;}
            if(ix!= NULL)
            {   free(ix); ix = NULL;}
            if(o2list!= NULL)
            {   free(o2list); o2list = NULL;}
            if(dyx!= NULL)
            {   free(dyx);dyx = NULL;}

        }      

    }
    
    /*new pid in the inference stage*/
    newid = 0;
    // allocate memory for elements
    output[0].pinstances[0].elements = (double*)malloc(ncount*dimsOutput*sizeof(double));
    output[0].lambda[0].elements  =  (struct cellData*)malloc(ncount*sizeof(struct cellData));
    col = dimsOutput;
    // copy results from the point
    for(iinst =0; iinst < ninst; iinst++)
        if(indicator[iinst]){ //update a match
                
                output[0].pinstances[0].elements[px2(newid, 3, ncount, col)] = point[iinst].score; //matched score
                output[0].pinstances[0].elements[px2(newid, 0, ncount, col)] = point[iinst].compid;// matched top level partid
                output[0].pinstances[0].elements[px2(newid, 1, ncount, col)] = pinstances[px2(iinst, 1, ninst, col)];//y
                output[0].pinstances[0].elements[px2(newid, 2, ncount, col)] = pinstances[px2(iinst, 2, ninst, col)];//x
                //allocate memory for its children
                output[0].lambda[0].elements[newid].a = (double*)malloc(pDim*point[iinst].length*sizeof(double));// y and x: two dimensions
                for(j=0; j< pDim*point[iinst].length; j++)
                    output[0].lambda[0].elements[newid].a[j] =  point[iinst].a[j];
                output[0].lambda[0].elements[newid].length = point[iinst].length;
                newid = newid+1;
                
    }
    
    output[0].pinstances[0].row = ncount;
    output[0].pinstances[0].col = dimsOutput;
    output[0].lambda[0].numel = ncount;
    if(indicator!=NULL)
    {  
        free(indicator); 
        indicator = NULL;
    }
    if(point!=NULL)
    {    free(point); point = NULL; }
        
    if(pdist!=NULL)
    {
        free(pdist);pdist=NULL;
    }
    if(pinstances!=NULL)
    {   free(pinstances);
        pinstances = NULL;
    }
    
    //return output;
}


/*fast inference on same image (different scales) with considering rotation for match*/
void comptopdown_fast(Data *output, unsigned char* morients, int imh, int imw, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, double rotation)
{
    
     //Data *output = (Data*)malloc(sizeof(Data));
     int i,j,k;
     int id;
     int ninst, col, iinst, newid, ncount;
     int r, c, o1, o2;
     struct array1d idx, idx1;
     int pDim = 4;
     ninst = 0;
     /*get the number of sampled points*/
     for(i=0; i< imh; i++)
         for(j=0;j<imw;j++)
             if(morients[px2(i,j,imh,imw)] >0)
             ninst++;
     /*construct pinstances*/ 
     double *pinstances = (double*)malloc(pDim*ninst*(sizeof(double)));
     id = 0; 
     for(i=0; i< imh; i++)
         for(j=0;j<imw;j++){
             if(morients[px2(i,j,imh,imw)] >0)
             {
             pinstances[px2(id, 0, ninst, pDim)] = morients[px2(i,j, imh, imw)];
             pinstances[px2(id, 1, ninst, pDim)] = i+1;//y coordinate
             pinstances[px2(id, 2, ninst, pDim)] = j+1;//x coordinate
             pinstances[px2(id, 3, ninst, pDim)] = 1;//score or flag
             id = id+1;}
     }

    /*define the output and allocate space*/
    int dimsOutput = 4;
    // allocate memory
    output[0].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
    output[0].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

    if(ninst==0)
    {
        output[0].pinstances[0].row = 0;
        output[0].pinstances[0].col = 0;
        output[0].lambda[0].numel =0;
        return;
    }
    
    int radius = model[numlayer-1].radius;//get the radius from the model
    int pre_radius = 0;
    bool *indicator = (bool*)malloc(ninst*sizeof(bool));
    for(iinst =0; iinst < ninst; iinst++)
        indicator[iinst] = false;
    struct Node *point = (struct Node*)malloc(ninst*sizeof(struct Node));
    /*counting the matched shaped*/
    ncount = 0;
    //double *pdist = funDist(pinstances, ninst, pDim);
    int *pdist = (int*)malloc(ninst*ninst*sizeof(int));
    funDist(pinstances, ninst, col, pdist); //compute the distance
    //#pragma omp parallel for 
    for(iinst = 0; iinst< ninst; iinst++) //match model by window scanning
    {


        pre_radius = 4;
        if(numlayer>1)
            pre_radius = model[numlayer-2].radius - radjust[numlayer-2];
        idx = searchneighs(iinst, pdist, ninst, radius, pre_radius, numlayer);
        if(idx.length>0){
            bool *visited = (bool*)malloc(ninst*sizeof(bool));
            initialize(visited, ninst);
            /*bool found = topdownmatch_fast(pinstances, ninst, pDim, pdist, iinst, idx, visited, model, numlayer, numOrient, 
                        numsub, radjust, discount,rotation, &point[iinst]);*/
            bool found=false;
            found = exactmatch(pinstances, ninst, pDim, pdist, iinst, idx, visited, model, numlayer, numOrient, numsub, 
                    radjust, discount,rotation, &point[iinst]);
            indicator[iinst] = found;
            if(found)     /*counting in the inference stage*/
                ncount++;
            free(visited);
        }      

    }
    
    /*new pid in the inference stage*/
    newid = 0;
    // allocate memory for elements
    output[0].pinstances[0].elements = (double*)malloc(ncount*dimsOutput*sizeof(double));
    output[0].lambda[0].elements  =  (struct cellData*)malloc(ncount*sizeof(struct cellData));
    col = dimsOutput;
    // copy results from the point
    for(iinst =0; iinst < ninst; iinst++)
        if(indicator[iinst]){ //update a match
                
                output[0].pinstances[0].elements[px2(newid, 3, ncount, col)] = point[iinst].score; //matched score
                output[0].pinstances[0].elements[px2(newid, 0, ncount, col)] = point[iinst].compid;// matched top level partid
                output[0].pinstances[0].elements[px2(newid, 1, ncount, col)] = pinstances[px2(iinst, 1, ninst, col)];//y
                output[0].pinstances[0].elements[px2(newid, 2, ncount, col)] = pinstances[px2(iinst, 2, ninst, col)];//x
                //allocate memory for its children
                output[0].lambda[0].elements[newid].a = (double*)malloc(pDim*point[iinst].length*sizeof(double));// y and x: two dimensions
                for(j=0; j< pDim*point[iinst].length; j++)
                    output[0].lambda[0].elements[newid].a[j] =  point[iinst].a[j];
                output[0].lambda[0].elements[newid].length = point[iinst].length;
                newid = newid+1;
                
    }
    
    output[0].pinstances[0].row = ncount;
    output[0].pinstances[0].col = dimsOutput;
    output[0].lambda[0].numel = ncount;
    if(indicator!=NULL)
    {  
        free(indicator); 
        indicator = NULL;
    }
    if(point!=NULL)
    {    free(point); point = NULL; }
    if(pdist!=NULL)
    {
        free(pdist);pdist=NULL;
    }
    if(pinstances!=NULL)
    {   free(pinstances);
        pinstances = NULL;
    }
    
    //return output;
}

/* inference on same image (different scales) with considering rotation for match*/
void comptopdown_full(Data *output, unsigned char* morients, int imh, int imw, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, double rotation)
{
    
     //Data *output = (Data*)malloc(sizeof(Data));
     int i,j,k;
     int id;
     int ninst, col, iinst, newid, ncount;
     int r, c, o1, o2;
     struct array1d idx, idx1;
     int pDim = 4;
     ninst = 0;
     /*get the number of sampled points*/
     for(i=0; i< imh; i++)
         for(j=0;j<imw;j++)
             if(morients[px2(i,j,imh,imw)] >0)
             ninst++;
     /*construct pinstances*/ 
     double *pinstances = (double*)malloc(pDim*ninst*(sizeof(double)));
     id = 0; 
     for(i=0; i< imh; i++)
         for(j=0;j<imw;j++){
             if(morients[px2(i,j,imh,imw)] >0)
             {
             pinstances[px2(id, 0, ninst, pDim)] = morients[px2(i,j, imh, imw)];
             pinstances[px2(id, 1, ninst, pDim)] = i+1;//y coordinate
             pinstances[px2(id, 2, ninst, pDim)] = j+1;//x coordinate
             pinstances[px2(id, 3, ninst, pDim)] = 1;//score or flag
             id = id+1;}
     }

    /*define the output and allocate space*/
    int dimsOutput = 4;
    // allocate memory
    output[0].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
    output[0].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

    if(ninst==0)
    {
        output[0].pinstances[0].row = 0;
        output[0].pinstances[0].col = 0;
        output[0].lambda[0].numel =0;
        return;
    }
    
    int radius = model[numlayer-1].radius;//get the radius from the model
    int pre_radius = 0;
    bool *indicator = (bool*)malloc(ninst*sizeof(bool));
    for(iinst =0; iinst < ninst; iinst++)
        indicator[iinst] = false;
    struct Node *point = (struct Node*)malloc(ninst*sizeof(struct Node));
    /*counting the matched shaped*/
    ncount = 0;
    //double *pdist = funDist(pinstances, ninst, pDim);
    int *pdist = (int*)malloc(ninst*ninst*sizeof(int));
    funDist(pinstances, ninst, col, pdist);
    //#pragma omp parallel for 
    for(iinst = 0; iinst< ninst; iinst++) //match model by window scanning
    {


        pre_radius = 4;
        if(numlayer>1)
            pre_radius = model[numlayer-2].radius - radjust[numlayer-2];
        idx = searchneighs(iinst, pdist, ninst, radius, pre_radius, numlayer);
        if(idx.length>0){
            bool *visited = (bool*)malloc(ninst*sizeof(bool));
            initialize(visited, ninst);
            bool found = topdownmatch_full(pinstances, ninst, pDim, pdist, iinst, idx, visited, model, numlayer, numOrient, 
                        numsub, radjust, discount,rotation, &point[iinst]);
            /*bool found=false;
            found = exactmatch(pinstances, ninst, pDim, pdist, iinst, idx, visited, model, numlayer, numOrient, numsub, 
                    radjust, discount,rotation, &point[iinst]);*/
            indicator[iinst] = found;
            if(found)     /*counting in the inference stage*/
                ncount++;
            free(visited);
        }      

    }
    
    /*new pid in the inference stage*/
    newid = 0;
    // allocate memory for elements
    output[0].pinstances[0].elements = (double*)malloc(ncount*dimsOutput*sizeof(double));
    output[0].lambda[0].elements  =  (struct cellData*)malloc(ncount*sizeof(struct cellData));
    col = dimsOutput;
    // copy results from the point
    for(iinst =0; iinst < ninst; iinst++)
        if(indicator[iinst]){ //update a match
                
                output[0].pinstances[0].elements[px2(newid, 3, ncount, col)] = point[iinst].score; //matched score
                output[0].pinstances[0].elements[px2(newid, 0, ncount, col)] = point[iinst].compid;// matched top level partid
                output[0].pinstances[0].elements[px2(newid, 1, ncount, col)] = pinstances[px2(iinst, 1, ninst, col)];//y
                output[0].pinstances[0].elements[px2(newid, 2, ncount, col)] = pinstances[px2(iinst, 2, ninst, col)];//x
                //allocate memory for its children
                output[0].lambda[0].elements[newid].a = (double*)malloc(pDim*point[iinst].length*sizeof(double));// y and x: two dimensions
                for(j=0; j< pDim*point[iinst].length; j++)
                    output[0].lambda[0].elements[newid].a[j] =  point[iinst].a[j];
                output[0].lambda[0].elements[newid].length = point[iinst].length;
                newid = newid+1;
                
    }
    
    output[0].pinstances[0].row = ncount;
    output[0].pinstances[0].col = dimsOutput;
    output[0].lambda[0].numel = ncount;
    if(indicator!=NULL)
    {  
        free(indicator); 
        indicator = NULL;
    }
    if(point!=NULL)
    {    free(point); point = NULL; }
        
    if(pdist!=NULL)
    {
        free(pdist);pdist=NULL;
    }
    if(pinstances!=NULL)
    {   free(pinstances);
        pinstances = NULL;
    }
    
    //return output;
}



/*group/batch inference with top down*/
Data* comptopdown_ex(struct Data *examples, int ndata, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust)
{
    
     Data *output = (Data*)malloc(sizeof(Data)*ndata);
     int i,j,k;
     int imidx,imh,imw,id;
     int ninst, col, iinst, newid, ncount;
     int r, c, o1, o2;
     struct array1d idx, idx1;
     int pDim = 4;
    
     for(imidx=0; imidx< ndata; imidx++) 
     {
        output[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
        output[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

     }
     int **pdist = (int**)malloc(ndata*sizeof(int*));
     /*inference first*/ 
     #pragma omp parallel for private(imidx)
     for(imidx=0; imidx< ndata; imidx++)
     {
         double *pinstances, *morients;
         imh = examples[imidx].imh;
         imw = examples[imidx].imw;
         ninst = 0;
         morients = examples[imidx].morients;
         if(examples[imidx].pinstances!=NULL && examples[imidx].pinstances[0].elements!=NULL){
             pinstances = examples[imidx].pinstances[0].elements;
             ninst = examples[imidx].pinstances[0].row;
             pDim = examples[imidx].pinstances[0].col;
             col = examples[imidx].pinstances[0].col;
         }
         else
         {    
             /*get the number of sampled points*/
             #pragma omp parallel for private(i,j)
             for(i=0; i< imh; i++)
                 for(j=0;j<imw;j++)
                     if(morients[px2(i,j,imh,imw)] >0)
                        ninst++;
             /*construct pinstances*/ 
             pinstances = (double*)malloc(pDim*ninst*(sizeof(double)));
             id = 0; 
             #pragma omp parallel for private(i,j)
             for(i=0; i< imh; i++)
                 for(j=0;j<imw;j++){
                     if(morients[px2(i,j,imh,imw)] >0)
                     {
                     pinstances[px2(id, 0, ninst, pDim)] = morients[px2(i,j, imh, imw)];
                     pinstances[px2(id, 1, ninst, pDim)] = i+1;//y coordinate
                     pinstances[px2(id, 2, ninst, pDim)] = j+1;//x coordinate
                     pinstances[px2(id, 3, ninst, pDim)] = 1;//score
                     id = id+1;}
             }
         }
        /*define the output and allocate space*/
        int dimsOutput = 4;
        // allocate memory
        //output[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
        //output[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

        if(ninst==0)
        {
            output[imidx].pinstances[0].row = 0;
            output[imidx].pinstances[0].col = 0;
            output[imidx].lambda[0].numel =0;
            pdist[imidx] = NULL;
            //return output;
            continue;
        }

        int radius = model[numlayer-1].radius;//get the radius from the model
        int pre_radius = 0;
        bool *indicator = (bool*)malloc(ninst*sizeof(bool));//initialize it
        for(iinst =0; iinst < ninst; iinst++)
             indicator[iinst] = false;
        struct Node *point = (struct Node*)malloc(ninst*sizeof(struct Node));
        // compute the distance
        //pdist[imidx] = funDist(pinstances, ninst, pDim);
        pdist[imidx] = (int*)malloc(ninst*ninst*sizeof(int));
        funDist(pinstances, ninst, col, pdist[imidx]);
        /*counting the matched shaped*/
        ncount = 0;
        #pragma omp parallel for private(iinst)
        for(iinst = 0; iinst< ninst; iinst++) //match model by window scanning
        {
            r = pinstances[px2(iinst, 1, ninst, pDim)];
            c = pinstances[px2(iinst, 2, ninst, pDim)];
            o1 = pinstances[px2(iinst, 0, ninst, pDim)];
            // find the neighbors 
            // idx = searchneighs(r, c, pinstances, ninst, pDim, radius, numlayer);
            pre_radius = 2;
            if(numlayer>1)
                pre_radius = model[numlayer-2].radius;
            // idx1 = searchneighs2(r, c, pinstances, ninst, pDim, pre_radius, numlayer-1, radjust);
            // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
            // idx = setdiff(idx,idx1);
            idx = searchneighs(iinst, pdist[imidx], ninst, radius, pre_radius, numlayer);
            if(idx.length>0)
            {
                int *y = (int*)malloc(sizeof(int)*idx.length); 
                //int x[idx.length];
                int *x = (int*)malloc(sizeof(int)*idx.length); 
                int *iy = (int*)malloc(sizeof(int)*idx.length); 
                int *ix = (int*)malloc(sizeof(int)*idx.length); 
                int *o2list = (int*)malloc(sizeof(int)*idx.length); 
                int *dyx = (int*)malloc(2*idx.length*sizeof(int));
                for(i=0;i< idx.length; i++)
                {

                    y[i] = (int)pinstances[px2(idx.a[i], 1, ninst, pDim)];//pinstances[idx.a[i]][2];
                    x[i] = (int)pinstances[px2(idx.a[i], 2, ninst, pDim)];//pinstances[idx.a[i]][3]; 
                    o2list[i] = (int)pinstances[px2(idx.a[i], 0, ninst, pDim)];//pinstances[idx.a[i]][1];
                    iy[i] = y[i] - r + radius;
                    ix[i] = x[i] -c + radius;
                    dyx[2*i]= iy[i]- radius;
                    dyx[2*i+1] = ix[i]-radius; 

                }

                int numdyx = idx.length;

                bool found = topdownmatch(pinstances, ninst, pDim, o1, r, c, dyx, o2list, numdyx, model, numlayer, numOrient, 
                        numsub, radjust, discount, &point[iinst]);


                indicator[iinst] = found;
                if(found)     /*counting in the inference stage*/
                    ncount++;
                if(y!=NULL)
                {   free(y); y=NULL;}
                if(x!=NULL)
                {   free(x); x= NULL;}
                if(iy!=NULL)
                {   free(iy); iy= NULL;}
                if(ix!= NULL)
                {   free(ix); ix = NULL;}
                if(o2list!= NULL)
                {   free(o2list); o2list = NULL;}
                if(dyx!= NULL)
                {   free(dyx);dyx = NULL;}

            }      

        }
        //printf("the number of instances: %d, detected: %d\n", ninst, ncount);
        /*new pid in the inference stage*/
        newid = 0;
        // allocate memory for elements
        output[imidx].pinstances[0].elements = (double*)malloc(ncount*dimsOutput*sizeof(double));
        output[imidx].lambda[0].elements  =  (struct cellData*)malloc(ncount*sizeof(struct cellData));
        // copy results from the point
        for(iinst =0; iinst < ninst; iinst++)
            if(indicator[iinst]){ //update a match

                    output[imidx].pinstances[0].elements[px2(newid, 3, ncount, pDim)] = point[iinst].score; //matched score
                    output[imidx].pinstances[0].elements[px2(newid, 0, ncount, pDim)] = point[iinst].compid;// matched top level partid
                    output[imidx].pinstances[0].elements[px2(newid, 1, ncount, pDim)] = pinstances[px2(iinst, 1, ninst, pDim)];//y
                    output[imidx].pinstances[0].elements[px2(newid, 2, ncount, pDim)] = pinstances[px2(iinst, 2, ninst, pDim)];//x
                    //allocate memory for its children
                    output[imidx].lambda[0].elements[newid].a = (double*)malloc(2*point[iinst].length*sizeof(double));// y and x: two dimensions
                    for(j=0; j< 2*point[iinst].length; j++)
                        output[imidx].lambda[0].elements[newid].a[j] =  point[iinst].a[j];
                    output[imidx].lambda[0].elements[newid].length = point[iinst].length;
                    newid = newid+1;

        }

        output[imidx].pinstances[0].row = ncount;
        output[imidx].pinstances[0].col = dimsOutput;
        output[imidx].lambda[0].numel = ncount;
        if(indicator!=NULL)
        {  
            free(indicator); 
            indicator = NULL;
        }
        if(point!=NULL)
        {    free(point); point = NULL; }

        if(pinstances!=NULL && (examples[imidx].pinstances==NULL || examples[imidx].pinstances[0].elements==NULL))
        {   free(pinstances);
            pinstances = NULL;
        }
    }
    //free the allocation for distance computation
    for(imidx=0; imidx<ndata; imidx++)
    {
        if(pdist[imidx]!=NULL)
        {   free(pdist[imidx]);
            pdist[imidx] = NULL;
        }
    }
    free(pdist);
     
    return output;
}



/*group/batch inference with top down*/
Data* comptopdown_ex(struct Data *examples, int ndata, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, int rotation)
{
    
     Data *output = (Data*)malloc(sizeof(Data)*ndata);
     int i,j,k;
     int imidx,imh,imw,id;
     int ninst, col, iinst, newid, ncount;
     int r, c, o1, o2;
     struct array1d idx, idx1;
     int pDim = 4;
    
     for(imidx=0; imidx< ndata; imidx++) 
     {
        output[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
        output[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

     }
     //allocate space for distance computing
     //double **pdist = (double**)malloc(ndata*sizeof(double*));
     int **pdist = (int**)malloc(ndata*sizeof(int*));
     /*inference first*/ 
     for(imidx=0; imidx< ndata; imidx++)
     {
         double *pinstances, *morients;
         imh = examples[imidx].imh;
         imw = examples[imidx].imw;
         ninst = 0;
         morients = examples[imidx].morients;
         if(examples[imidx].pinstances!=NULL && examples[imidx].pinstances[0].elements!=NULL){
             pinstances = examples[imidx].pinstances[0].elements;
             ninst = examples[imidx].pinstances[0].row;
             pDim = examples[imidx].pinstances[0].col;
             col = examples[imidx].pinstances[0].col;
         }
         else
         {    
             /*get the number of sampled points*/
             for(i=0; i< imh; i++)
                 for(j=0;j<imw;j++)
                     if(morients[px2(i,j,imh,imw)] >0)
                        ninst++;
             /*construct pinstances*/ 
             pinstances = (double*)malloc(pDim*ninst*(sizeof(double)));
             id = 0; 
             for(i=0; i< imh; i++)
                 for(j=0;j<imw;j++){
                     if(morients[px2(i,j,imh,imw)] >0)
                     {
                     pinstances[px2(id, 0, ninst, pDim)] = morients[px2(i,j, imh, imw)];
                     pinstances[px2(id, 1, ninst, pDim)] = i+1;//y coordinate
                     pinstances[px2(id, 2, ninst, pDim)] = j+1;//x coordinate
                     pinstances[px2(id, 3, ninst, pDim)] = 1;//score
                     id = id+1;}
             }
         }
        /*define the output and allocate space*/
        int dimsOutput = 4;
        // allocate memory
        //output[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
        //output[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

        if(ninst==0)
        {
            output[imidx].pinstances[0].row = 0;
            output[imidx].pinstances[0].col = 0;
            output[imidx].lambda[0].numel =0;
            pdist[imidx] = NULL;
            //return output;
            continue;
        }

        int radius = model[numlayer-1].radius;//get the radius from the model
        int pre_radius = 0;
        bool *indicator = (bool*)malloc(ninst*sizeof(bool));//initialize it
        for(iinst =0; iinst < ninst; iinst++)
             indicator[iinst] = false;
        struct Node *point = (struct Node*)malloc(ninst*sizeof(struct Node));
        //pdist[imidx] = funDist(pinstances, ninst, pDim);
        pdist[imidx] = (int*)malloc(ninst*ninst*sizeof(int));
        funDist(pinstances, ninst, col, pdist[imidx]);
        /*counting the matched shaped*/
        ncount = 0;
        for(iinst = 0; iinst< ninst; iinst++) //match model by window scanning
        {
            r = pinstances[px2(iinst, 1, ninst, pDim)];
            c = pinstances[px2(iinst, 2, ninst, pDim)];
            o1 = pinstances[px2(iinst, 0, ninst, pDim)];
            // find the neighbors 
            //idx = searchneighs(r, c, pinstances, ninst, pDim, radius, numlayer);
            pre_radius = 2;
            if(numlayer>1)
                pre_radius = model[numlayer-2].radius;
            // idx1 = searchneighs2(r, c, pinstances, ninst, pDim, pre_radius, numlayer-1, radjust);
            // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
            // idx = setdiff(idx,idx1);
            idx = searchneighs(iinst, pdist[imidx], ninst, radius, pre_radius, numlayer);
            if(idx.length>0)
            {
                int *y = (int*)malloc(sizeof(int)*idx.length); 
                //int x[idx.length];
                int *x = (int*)malloc(sizeof(int)*idx.length); 
                int *iy = (int*)malloc(sizeof(int)*idx.length); 
                int *ix = (int*)malloc(sizeof(int)*idx.length); 
                int *o2list = (int*)malloc(sizeof(int)*idx.length); 
                int *dyx = (int*)malloc(2*idx.length*sizeof(int));
                for(i=0;i< idx.length; i++)
                {

                    y[i] = (int)pinstances[px2(idx.a[i], 1, ninst, pDim)];//pinstances[idx.a[i]][2];
                    x[i] = (int)pinstances[px2(idx.a[i], 2, ninst, pDim)];//pinstances[idx.a[i]][3]; 
                    o2list[i] = (int)pinstances[px2(idx.a[i], 0, ninst, pDim)];//pinstances[idx.a[i]][1];
                    iy[i] = y[i] - r + radius;
                    ix[i] = x[i] -c + radius;
                    dyx[2*i]= iy[i]- radius;
                    dyx[2*i+1] = ix[i]-radius; 

                }

                int numdyx = idx.length;

                bool found = topdownmatch(pinstances, ninst, pDim, o1, r, c, dyx, o2list, numdyx, model, numlayer, numOrient, 
                        numsub, radjust, discount, &point[iinst]);


                indicator[iinst] = found;
                if(found)     /*counting in the inference stage*/
                    ncount++;
                if(y!=NULL)
                {   free(y); y=NULL;}
                if(x!=NULL)
                {   free(x); x= NULL;}
                if(iy!=NULL)
                {   free(iy); iy= NULL;}
                if(ix!= NULL)
                {   free(ix); ix = NULL;}
                if(o2list!= NULL)
                {   free(o2list); o2list = NULL;}
                if(dyx!= NULL)
                {   free(dyx);dyx = NULL;}

            }      

        }
        //printf("the number of instances: %d, detected: %d\n", ninst, ncount);
        /*new pid in the inference stage*/
        newid = 0;
        // allocate memory for elements
        output[imidx].pinstances[0].elements = (double*)malloc(ncount*dimsOutput*sizeof(double));
        output[imidx].lambda[0].elements  =  (struct cellData*)malloc(ncount*sizeof(struct cellData));
        // copy results from the point
        for(iinst =0; iinst < ninst; iinst++)
            if(indicator[iinst]){ //update a match

                    output[imidx].pinstances[0].elements[px2(newid, 3, ncount, pDim)] = point[iinst].score; //matched score
                    output[imidx].pinstances[0].elements[px2(newid, 0, ncount, pDim)] = point[iinst].compid;// matched top level partid
                    output[imidx].pinstances[0].elements[px2(newid, 1, ncount, pDim)] = pinstances[px2(iinst, 1, ninst, pDim)];//y
                    output[imidx].pinstances[0].elements[px2(newid, 2, ncount, pDim)] = pinstances[px2(iinst, 2, ninst, pDim)];//x
                    //allocate memory for its children
                    output[imidx].lambda[0].elements[newid].a = (double*)malloc(2*point[iinst].length*sizeof(double));// y and x: two dimensions
                    for(j=0; j< 2*point[iinst].length; j++)
                        output[imidx].lambda[0].elements[newid].a[j] =  point[iinst].a[j];
                    output[imidx].lambda[0].elements[newid].length = point[iinst].length;
                    newid = newid+1;

        }

        output[imidx].pinstances[0].row = ncount;
        output[imidx].pinstances[0].col = dimsOutput;
        output[imidx].lambda[0].numel = ncount;
        if(indicator!=NULL)
        {  
            free(indicator); 
            indicator = NULL;
        }
        if(point!=NULL)
        {    free(point); point = NULL; }

        if(pinstances!=NULL && (examples[imidx].pinstances==NULL || examples[imidx].pinstances[0].elements==NULL))
        {   free(pinstances);
            pinstances = NULL;
        }
    }
     
    //free the allocation for distance computation
    for(imidx=0; imidx<ndata; imidx++)
    {
        if(pdist[imidx]!=NULL)
        {   free(pdist[imidx]);
            pdist[imidx] = NULL;
        }
    }
    free(pdist);
     
    return output;
}


/*group/batch inference with top down*/
Data* comptopdown_ex2(struct Data *examples, int ndata, struct Parts *model, int numlayer, int numsub, 
        int numOrient, double discount, double *radjust, int rotation)
{
    
     Data *output = (Data*)malloc(sizeof(Data)*ndata);
     int i,j,k;
     int imidx,imh,imw,id;
     int ninst, col, iinst, newid, ncount;
     int r, c, o1, o2;
     struct array1d idx, idx1;
     int pDim = 4;
    
     for(imidx=0; imidx< ndata; imidx++) 
     {
        output[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
        output[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

     }
     //allocate space for distance computing
     int **pdist = (int**)malloc(ndata*sizeof(int*));
     /*inference first*/ 
     for(imidx=0; imidx< ndata; imidx++)
     {
         double *pinstances, *morients;
         imh = examples[imidx].imh;
         imw = examples[imidx].imw;
         ninst = 0;
         morients = examples[imidx].morients;
         if(examples[imidx].pinstances!=NULL && examples[imidx].pinstances[0].elements!=NULL){
             pinstances = examples[imidx].pinstances[0].elements;
             ninst = examples[imidx].pinstances[0].row;
             pDim = examples[imidx].pinstances[0].col;
             col = examples[imidx].pinstances[0].col;
         }
         else
         {    
             /*get the number of sampled points*/
             for(i=0; i< imh; i++)
                 for(j=0;j<imw;j++)
                     if(morients[px2(i,j,imh,imw)] >0)
                        ninst++;
             /*construct pinstances*/ 
             pinstances = (double*)malloc(pDim*ninst*(sizeof(double)));
             id = 0; 
             for(i=0; i< imh; i++)
                 for(j=0;j<imw;j++){
                     if(morients[px2(i,j,imh,imw)] >0)
                     {
                     pinstances[px2(id, 0, ninst, pDim)] = morients[px2(i,j, imh, imw)];
                     pinstances[px2(id, 1, ninst, pDim)] = i+1;//y coordinate
                     pinstances[px2(id, 2, ninst, pDim)] = j+1;//x coordinate
                     pinstances[px2(id, 3, ninst, pDim)] = 1;//score
                     id = id+1;}
             }
         }
        /*define the output and allocate space*/
        int dimsOutput = 4;
        // allocate memory
        //output[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
        //output[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

        if(ninst==0)
        {
            output[imidx].pinstances[0].row = 0;
            output[imidx].pinstances[0].col = 0;
            output[imidx].lambda[0].numel =0;
            pdist[imidx] = NULL;
            //return output;
            continue;
        }

        int radius = model[numlayer-1].radius;//get the radius from the model
        int pre_radius = 0;
        bool *indicator = (bool*)malloc(ninst*sizeof(bool));//initialize it
        for(iinst =0; iinst < ninst; iinst++)
             indicator[iinst] = false;
        struct Node *point = (struct Node*)malloc(ninst*sizeof(struct Node));
        //compute distances
        //pdist[imidx] = funDist(pinstances, ninst, pDim);
        pdist[imidx] = (int*)malloc(ninst*ninst*sizeof(int));
        funDist(pinstances, ninst, col, pdist[imidx]);
        /*counting the matched shaped*/
        ncount = 0;
        for(iinst = 0; iinst< ninst; iinst++) //match model by window scanning
        {
            r = pinstances[px2(iinst, 1, ninst, pDim)];
            c = pinstances[px2(iinst, 2, ninst, pDim)];
            o1 = pinstances[px2(iinst, 0, ninst, pDim)];
            // find the neighbors 
            //idx = searchneighs(r, c, pinstances, ninst, pDim, radius, numlayer);
            pre_radius = 2;
            if(numlayer>1)
                pre_radius = model[numlayer-2].radius;
            
            idx = searchneighs(iinst, pdist[imidx], ninst, radius, pre_radius, numlayer);
            if(idx.length>0)
            {   
                bool *visited = (bool*)malloc(ninst*sizeof(bool));
                initialize(visited, ninst);
                bool found = topdownmatch_fast(pinstances, ninst, pDim, pdist[imidx], iinst, idx, visited, model, numlayer, numOrient, 
                        numsub, radjust, discount,rotation, &point[iinst]);

                indicator[iinst] = found;
                if(found)     /*counting in the inference stage*/
                    ncount++;
                free(visited);

            }      

        }
        //printf("the number of instances: %d, detected: %d\n", ninst, ncount);
        /*new pid in the inference stage*/
        newid = 0;
        // allocate memory for elements
        output[imidx].pinstances[0].elements = (double*)malloc(ncount*dimsOutput*sizeof(double));
        output[imidx].lambda[0].elements  =  (struct cellData*)malloc(ncount*sizeof(struct cellData));
        // copy results from the point
        for(iinst =0; iinst < ninst; iinst++)
            if(indicator[iinst]){ //update a match

                    output[imidx].pinstances[0].elements[px2(newid, 3, ncount, pDim)] = point[iinst].score; //matched score
                    output[imidx].pinstances[0].elements[px2(newid, 0, ncount, pDim)] = point[iinst].compid;// matched top level partid
                    output[imidx].pinstances[0].elements[px2(newid, 1, ncount, pDim)] = pinstances[px2(iinst, 1, ninst, pDim)];//y
                    output[imidx].pinstances[0].elements[px2(newid, 2, ncount, pDim)] = pinstances[px2(iinst, 2, ninst, pDim)];//x
                    //allocate memory for its children
                    output[imidx].lambda[0].elements[newid].a = (double*)malloc(2*point[iinst].length*sizeof(double));// y and x: two dimensions
                    for(j=0; j< point[iinst].length; j++){
                        output[imidx].lambda[0].elements[newid].a[2*j+0] =  point[iinst].a[j*pDim+0];
                        output[imidx].lambda[0].elements[newid].a[2*j+1] =  point[iinst].a[j*pDim+1];
                    }
                    output[imidx].lambda[0].elements[newid].length = point[iinst].length;
                    newid = newid+1;

        }

        output[imidx].pinstances[0].row = ncount;
        output[imidx].pinstances[0].col = dimsOutput;
        output[imidx].lambda[0].numel = ncount;
        if(indicator!=NULL)
        {  
            free(indicator); 
            indicator = NULL;
        }
        if(point!=NULL)
        {    free(point); point = NULL; }

        if(pinstances!=NULL && (examples[imidx].pinstances==NULL || examples[imidx].pinstances[0].elements==NULL))
        {   free(pinstances);
            pinstances = NULL;
        }
    }
     
    //free the allocation for distance computation
    for(imidx=0; imidx<ndata; imidx++)
    {
        if(pdist[imidx]!=NULL)
        {   free(pdist[imidx]);
            pdist[imidx] = NULL;
        }
    }
    free(pdist);
     
    return output;
}


/*
struct array1d matchmid(int o1, struct Parts *model,int numlayer, int numsub)
{
    struct array1d reidx;
    int count = 0;
    int ntypes = model[numlayer].ptypes_[numsub-1].ntypes;
    if(ntypes != model[numlayer].ntypes_[numsub-1])
    {    
        printf("the number of parts is not consistent!\n");
        exit(0);
    }
    for(int i=0; i<ntypes; i++)
    {
        if(model[numlayer].ptypes_[numsub-1].dyx[i].a[1] == o1)
        {
            //save pid here//
            reidx.a[count] = model[numlayer].ptypes_[numsub-1].dyx[i].a[0];
            count = count + 1;
        }
    
    }
    reidx.length = count;
    return reidx;
    
}*/

/* -----------------------------------------------------------------------------
 *                        define some functions 
 * -----------------------------------------------------------------------------
 */      


int getnumchild(int centerid, struct Parts *model, int numlayer, int numsub)
{
    int numchild = 0;
    int allchilds = 0;
    int anchorid, childid;

    centerid = centerid -1; //index from 0, so minus 1
    int lid= numlayer;//-1;
    //assert(lid>1);
    if(lid>1){
            
        numchild = model[lid-1].ptypes_[numsub-1].dyx[centerid].a[2];
        //add itself
        allchilds = allchilds + numchild;

        //consider the center id
        anchorid = model[lid-1].ptypes_[numsub-1].dyx[centerid].a[1];
        //centerid = centerid -1;
        //anchor id for recursive
        allchilds = allchilds + getnumchild(anchorid, model, lid-1, numsub);

        //add its children
        for(int i= 0; i< numchild; i++){       
            childid = model[lid-1].ptypes_[numsub-1].dyx[centerid].a[3+5*i];
            //childid = childid-1; // index starts from zeros
            //recursive to previous layer
            allchilds = allchilds + getnumchild(childid, model, lid-1, numsub);
            
        }

    }
    else
    {    return 0;
    }
    
    //centerid = centerid + 1;//restore it back to 1: numOrient
    return allchilds;
}


int id2orient(int centerid, struct Parts *model, int numlayer, int numsub)
{
    centerid = centerid -1; //index from 0, so minus 1
    for(int lid= numlayer-1; lid>1; lid--)
    {
        centerid = model[lid-1].ptypes_[numsub-1].dyx[centerid].a[1];
        centerid = centerid -1;
    }
    centerid = centerid + 1;//restore it back to 1: numOrient
    return centerid;
}


//initialize img data
void intializeimg(int *img, int imh, int imw)
{
int i,j;
for(i=0; i< imh; i++)
for(j=0; j< imw; j++)
{   
   
   img[px2(i,j,imh,imw)] = 0;

}
}


//copy filter data to img data from position y and x
void copyimg(int *img, int imh, int imw, int *filter, int h,int w, int y, int x)
{
int i,ii;
int j,jj;
int miny = MAX(y,1);
int minx = MAX(x,1);
int maxy = MIN(imh, y+h);
int maxx = MIN(imw, x+w);
for(i=miny; i< maxy; i++)
for(j=minx; j< maxx; j++)
{   
   ii = i-miny; jj = j-minx;
   img[px2(i,j,imh,imw)] = filter[px2(ii,jj,h,w)];

}
}

//copy filter data to img data from position y and x
void addimg(int *img, int imh, int imw, int *filter, int h,int w, int y, int x)
{
int i,ii;
int j,jj;
int miny = MAX(y,1);
int minx = MAX(x,1);
int maxy = MIN(imh, y+h);
int maxx = MIN(imw, x+w);
for(i=miny; i< maxy; i++)
    for(j=minx; j< maxx; j++)
    {   
       ii = i-miny; jj = j-minx;
       img[px2(i,j,imh,imw)] =img[px2(i,j,imh,imw)] + filter[px2(ii,jj,h,w)];

    }
}


double* drawsubparts(int *img, int imh, int imw, int y, int x,   struct Parts *model,int centerid, int numsub, int numOrient, int numlayer, int **allFilter, int h, int w)
{

    int anchorid, childid, numchild;
    int dy, dx, cy,cx;

    anchorid = model[numlayer-1].ptypes_[numsub-1].dyx[centerid].a[1];
    int ori1 = id2orient(anchorid, model, numlayer, numsub);
    // add the filter to the reconstructed image
    addimg(img, imh, imw, allFilter[ori1], h, w, y, x);

    /* child of anchor id */
    if(numlayer>2)
        drawsubparts(img, imh, imw, y,x, model, anchorid, numsub, numOrient, numlayer-1, allFilter, w, h); 


    numchild = model[numlayer-1].ptypes_[numsub-1].dyx[centerid].a[2];
    // allocate space for its children
    int **subimg = (int**)malloc(numchild*sizeof(int*));
    for(int j=0; j<numchild;j++){
        subimg[j] = (int*)malloc(imh*imw*sizeof(int));
        intializeimg(subimg[j], imh,imw);
    }

    /* child of the current part id*/
    for(int j=0; j<numchild;j++){
         childid = model[numlayer-1].ptypes_[numsub-1].dyx[centerid].a[3+5*(j-1)];
         dy = model[numlayer-1].ptypes_[numsub-1].dyx[centerid].a[4+5*(j-1)];
         dx = model[numlayer-1].ptypes_[numsub-1].dyx[centerid].a[5+5*(j-1)];
         cy = y + dy; cy = floor(cy);
         cx = x + dx; cx = floor(cx);
         int ori2 = id2orient(childid, model, numlayer, numsub);
         addimg(subimg[j], imh, imw, allFilter[ori2], h, w, cy, cx);
         if(numlayer>2)
            drawsubparts(subimg[j], imh, imw, cy,cx, model, childid, numsub, numOrient, numlayer-1, allFilter, w, h); 

    }
    
    
    
    double *overlap = (double*)malloc(numchild*sizeof(double));
    // dilate the image
    for(int i=0;i<2;i++)
    {
        imdilate(img, imh, imw);
        for(int j=0; j<numchild;j++){
            imdilate(subimg[j], imh, imw);;
        }
    }
    //caculate whether the anchor part is overlap with its childrens
    for(int j=0; j<numchild;j++)
        overlap[j] = imoverlap(img, subimg[j], imh, imw);

    // get the final image
    for(int j=0; j<numchild;j++){
        addimg(img, imh,imw, subimg[j],imh,imw,0,0);
    }

    
    // delete space for its children
    for(int j=0; j<numchild;j++){
        free(subimg[j]);
    }
    free(subimg);
    return overlap;
}


struct pTypes* mappingParts(struct Parts *model, int numOrient, int numlayer, int numsub, struct array1d *tdlinks)
{
    int i,j;
    int pDim = 9;
    // struct pTypes *ptypes = model[numlayer-].ptypes[numsub-1].dyx;
    //double *ptypes = (double*)malloc(pDim*numparts*sizeof(double));
    struct pTypes *flatparts;
    if(numlayer<2)
    {
        flatparts = NULL;
        return flatparts;
    }
    flatparts = (pTypes*)malloc(sizeof(pTypes));
    int numparts = model[numlayer-1].ntypes_[numsub-1];
    int centerid, childid, numchild;
    for(i=0; i< numparts; i++)
    {    
        centerid = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[1];
        numchild = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[2];
        //copy parts information first
        flatparts[0].dyx[i].a[0] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[0];//pid
        //model[numlayer-1].ptypes_[numsub-1].dyx[i].a[0];//pid
        flatparts[0].dyx[i].a[1] = id2orient(centerid, model, numlayer, numsub);//anchor id
        flatparts[0].dyx[i].a[2] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[2];//#children
        // check its children
        for(j=0; j< numchild; j++){
            //get the flatten types
            childid = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[3+5*j];//childid
            flatparts[0].dyx[i].a[3+5*j] = id2orient(childid, model, numlayer, numsub);
            flatparts[0].dyx[i].a[4+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[4+5*j];;
            flatparts[0].dyx[i].a[5+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[5+5*j];
            flatparts[0].dyx[i].a[6+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[6+5*j];
            flatparts[0].dyx[i].a[7+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[7+5*j];
            // model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = ;      
            flatparts[0].dyx[i].a[8+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[8+5*j];  
        
        }
        flatparts[0].dyx[i].length = 3+ 5*numchild;
        flatparts[0].ntypes = numparts;
        
    }
    
    // get links from top down
    
    int *ncount = (int*)malloc(numOrient*sizeof(int));
    for(i=0;i<numOrient; i++)
        ncount[i] = 0;
    for(i=0; i< numparts; i++)
    {
        //len =length(flatparts[i].length);
        centerid = flatparts[0].dyx[i].a[1] -1; //start from 0
        tdlinks[centerid].a[ncount[centerid]] = flatparts[0].dyx[i].a[0] -1; //point to pid & start from zero
        ncount[centerid] = ncount[centerid]+1;
    }
    for(i=0;i<numOrient; i++)
        tdlinks[i].length = ncount[i];
    
    // unique each links
    for(i=0; i< numOrient; i++){
        tdlinks[i].length = ncount[i];
        //tdlinks[i] = unique(tdlinks[i], 1);
    }
    if(ncount!=NULL)
    {free(ncount);ncount =NULL;}
    return flatparts;
}



/*flatten parts with orientations*/
struct pTypes* mappingParts(struct Parts *model, int numOrient, int numlayer, int numsub, int rotation, struct array1d *tdlinks)
{
    int i,j,otemp;
    int pDim = 9;
    // struct pTypes *ptypes = model[numlayer-].ptypes[numsub-1].dyx;
    //double *ptypes = (double*)malloc(pDim*numparts*sizeof(double));
    struct pTypes *flatparts;
    if(numlayer<2)
    {
        flatparts = NULL;
        return flatparts;
    }
    flatparts = (pTypes*)malloc(sizeof(pTypes));  
    int numparts = model[numlayer-1].ntypes_[numsub-1];
    int centerid, childid, numchild;
    for(i=0; i< numparts; i++)
    {    
        centerid = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[1];
        numchild = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[2];
        //copy parts information first
        flatparts[0].dyx[i].a[0] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[0];//pid
        //model[numlayer-1].ptypes_[numsub-1].dyx[i].a[0];//pid
        otemp = id2orient(centerid, model, numlayer, numsub);//anchor id by flattern
        //add rotation
        otemp = otemp + rotation;
        if(otemp<1) otemp = otemp + numOrient;
        if(otemp>numOrient) otemp = otemp - numOrient;
        flatparts[0].dyx[i].a[1] = otemp; //new anchor id here
        flatparts[0].dyx[i].a[2] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[2];//#children
        // check its children
        for(j=0; j< numchild; j++){
            //get the flatten types
            childid = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[3+5*j];//childid
            otemp = id2orient(childid, model, numlayer, numsub);
            //add rotation
            otemp = otemp + rotation;
            if(otemp<1) otemp = otemp + numOrient;
            if(otemp>numOrient) otemp = otemp - numOrient;
            flatparts[0].dyx[i].a[3+5*j] = otemp;//new child id from flatten
            flatparts[0].dyx[i].a[4+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[4+5*j];;
            flatparts[0].dyx[i].a[5+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[5+5*j];
            flatparts[0].dyx[i].a[6+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[6+5*j];
            flatparts[0].dyx[i].a[7+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[7+5*j];
            // model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = ;      
            flatparts[0].dyx[i].a[8+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[8+5*j];  
        
        }
        flatparts[0].dyx[i].length = 3+ 5*numchild;
        flatparts[0].ntypes = numparts;
        
    }
    
    // get links from top down
    
    int *ncount = (int*)malloc(numOrient*sizeof(int));
    for(i=0;i<numOrient; i++)
        ncount[i] = 0;
    for(i=0; i< numparts; i++)
    {
        //len =length(flatparts[i].length);
        centerid = flatparts[0].dyx[i].a[1] -1; //start from 0
        tdlinks[centerid].a[ncount[centerid]] = flatparts[0].dyx[i].a[0] -1; //point to pid & start from zero
        ncount[centerid] = ncount[centerid]+1;
    }
    for(i=0;i<numOrient; i++)
        tdlinks[i].length = ncount[i];
    
    // unique each links
    for(i=0; i< numOrient; i++){
        tdlinks[i].length = ncount[i];
        //tdlinks[i] = unique(tdlinks[i], 1);
    }
    if(ncount!=NULL)
    {free(ncount);ncount =NULL;}
    return flatparts;
}


/*flatten all layer parts with orientations*/
struct fParts* flattenparts(struct Parts *model, int numOrient, int numlayer, int numsub, int rotation)
{
    int i,j,otemp;
    int pDim = 9;
    // struct pTypes *ptypes = model[numlayer-].ptypes[numsub-1].dyx;
    //double *ptypes = (double*)malloc(pDim*numparts*sizeof(double));
    struct fParts *flatparts;
    if(numlayer<2)
    {
        flatparts = NULL;
        return flatparts;
    }
    flatparts= (fParts*)malloc(numlayer*sizeof(fParts)); 
    //initialize the empty links;
    for(i=0;i<numlayer;i++)
        flatparts[i].tdlinks = NULL;
    
    // processing each layer
    while(numlayer>2)
    {
        int numparts = model[numlayer-1].ntypes_[numsub-1];
        int centerid, childid, numchild;
        for(i=0; i< numparts; i++)
        {    
            centerid = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[1];
            numchild = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[2];
            //copy parts information first
            flatparts[numlayer-1].dyx[i].a[0] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[0];//pid
            //model[numlayer-1].ptypes_[numsub-1].dyx[i].a[0];//pid
            otemp = id2orient(centerid, model, numlayer, numsub);//anchor id by flattern
            //add rotation
            otemp = otemp + rotation;
            if(otemp<1) otemp = otemp + numOrient;
            if(otemp>numOrient) otemp = otemp - numOrient;
            flatparts[numlayer-1].dyx[i].a[1] = otemp; //new anchor id here
            flatparts[numlayer-1].dyx[i].a[2] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[2];//#children
            // check its children
            for(j=0; j< numchild; j++){
                //get the flatten types
                childid = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[3+5*j];//childid
                otemp = id2orient(childid, model, numlayer, numsub);
                //add rotation
                otemp = otemp + rotation;
                if(otemp<1) otemp = otemp + numOrient;
                if(otemp>numOrient) otemp = otemp - numOrient;
                flatparts[numlayer-1].dyx[i].a[3+5*j] = otemp;//new child id from flatten
                flatparts[numlayer-1].dyx[i].a[4+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[4+5*j];;
                flatparts[numlayer-1].dyx[i].a[5+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[5+5*j];
                flatparts[numlayer-1].dyx[i].a[6+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[6+5*j];
                flatparts[numlayer-1].dyx[i].a[7+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[7+5*j];
                // model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = ;      
                flatparts[numlayer-1].dyx[i].a[8+5*j] = model[numlayer-1].ptypes_[numsub-1].dyx[i].a[8+5*j];  

            }
            flatparts[numlayer-1].dyx[i].length = 3+ 5*numchild;
            flatparts[numlayer-1].ntypes = numparts;

        }

        // get links from top down
        flatparts[numlayer-1].tdlinks = (struct array1d*)malloc(numOrient*sizeof(struct array1d));
        int *ncount = (int*)malloc(numOrient*sizeof(int));
        for(i=0;i<numOrient; i++)
            ncount[i] = 0;
        for(i=0; i< numparts; i++)
        {
            //len =length(flatparts[i].length);
            centerid = flatparts[numlayer-1].dyx[i].a[1] -1; //start from 0
            flatparts[numlayer-1].tdlinks[centerid].a[ncount[centerid]] = flatparts[numlayer-1].dyx[i].a[0] -1; //point to pid & start from zero
            ncount[centerid] = ncount[centerid]+1;
        }
        for(i=0;i<numOrient; i++)
            flatparts[numlayer-1].tdlinks[i].length = ncount[i];

        // unique each links
        for(i=0; i< numOrient; i++){
            flatparts[numlayer-1].tdlinks[i].length = ncount[i];
            //tdlinks[i] = unique(tdlinks[i], 1);
        }
        if(ncount!=NULL)
        {free(ncount);ncount =NULL;}
        
        numlayer = numlayer -1;
    }
    return flatparts;
}



struct array1d unique(struct array1d d, int dim)
{
    int i,j;
    assert(dim==1);
    int len = d.length;
    int *indicator = (int*)malloc(len*sizeof(int));
    for(i=0;i< d.length-1; i++)
        indicator[i] = 1;
    for(i=0;i< d.length-1; i++)
    {   
        if(indicator[i])
        for(j= i+1; j< len; j++)
            if(d.a[j]==d.a[i])
                indicator[j]= 0;
        
    }
    len = sum(indicator, d.length);
    struct array1d newd;
    j = 0;
    for(i=0;i< d.length-1; i++)
    {   
        if(indicator[i])
        {
            newd.a[j]= d.a[i];
            j = j+1;
        }
        
    }
    newd.length = j;
    if(indicator!=NULL)
    {free(indicator);indicator = NULL;}
    return newd;
}

void initialize(int *flag, int n){
    
	for(int i= 0; i<n; i++)
        flag[i] = 0;
   
}

void initialize(bool *flag, int n){
    
	for(int i= 0; i<n; i++)
        flag[i] = false;
   
}



int vmax(double *vector, int n, double *score){
  
    *score = -INF;
    int id=0;
    for(int i=0; i< n; i++)
       if(*score < vector[i])
       {    *score = vector[i];
            id = i;
       }
    return id;
}


struct Node cat(struct Node pchild, struct Node cchild)
{
    
    int len = MIN(pchild.length + cchild.length, MAX_LEN);
    int i;
    for(i=pchild.length; i< len; i++){
        pchild.a[2*i] = cchild.a[2*(i-pchild.length)];
        pchild.a[2*i + 1] = cchild.a[2*(i-pchild.length) +1];
    }
    pchild.length = len;
    return pchild;
}


struct Node cat(struct Node pchild, struct Node cchild, int pDim)
{
    
    int len = MIN(pchild.length + cchild.length, MAX_LEN);
    int i,j;
    for(i=pchild.length; i< len; i++){
        for(j=0; j<pDim;j++){
            pchild.a[pDim*i+j] = cchild.a[pDim*(i-pchild.length)+j];
        }
    }
    pchild.length = len;
    return pchild;
}

void cat(struct Node *pchild, struct Node *cchild, int pDim)
{
    
    int len = MIN(pchild->length + cchild->length, MAX_LEN);
    int i,j;
    for(i=pchild->length; i< len; i++){
        for(j=0; j<pDim;j++){
            pchild->a[pDim*i+j] = cchild->a[pDim*(i-pchild->length)+j];
        }
    }
    pchild->length = len;
    
}

/*
* ---------------------------------------------------------------
*            top down matching in an recursively way
* ---------------------------------------------------------------
/


/* WITH TOP DOWN MATCHING */
bool topdownmatch(double *pinstances, int ninst, int col, int o1, int r, int c, int *dyx, int *o2list,int numdyx,
        Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount, Node *point)
{   
    // initialize the value for return
    bool found, cfound;
    found= false; cfound = false;   

    struct array1d *tdlinks = (struct array1d*)malloc(numOrient*sizeof(struct array1d));
    struct pTypes *flatparts = mappingParts(model, numOrient, numlayer, numsub, tdlinks);
    if(flatparts==NULL)
    {   found=false;
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        return found;
    }
    //get the matching parts id
    struct array1d partids = tdlinks[o1-1];
    if(partids.length<1){
        found=false;
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        return found;
    } 
    
    // define variables
    int i,j,k;
    double epsilon = 1.5;
    double mean[2]; double var[2];
    int childid, numchild, anchorid;
    int cy, cx, co2, radius, pre_radius;//for its children
    struct array1d reidx, idx1;
    // initialize it
    int *indicator = (int*)malloc(partids.length*sizeof(int));
    double *tscores = (double*)malloc(partids.length*sizeof(double));
    int *childlens = (int*)malloc(partids.length*sizeof(int));
    for(i=0;i< partids.length; i++)
    {
        indicator[i] = 0;
        tscores[i] = 0;
        childlens[i] = 1;
    }
    //define the part nodes by the number of pids
    struct Node *partnodes = (struct Node*)malloc(partids.length*sizeof(struct Node));
    // check its neighborhood
    int tindex =-1;
    int nummatched = 0;
    double maxprob=0;
    for(i=0; i<partids.length; i++)
    { // processing each part type
        tindex = tindex + 1;
        partnodes[i].length = 0; //initialize it
        numchild = flatparts[0].dyx[partids.a[i]].a[2];
        childlens[i] = numchild;
        //indicator matching for child
        int* flag = (int*)malloc(numchild*sizeof(int));
        initialize(flag, numchild);
        //define the current node from this input point
        struct Node pnode;
        pnode.a[0] = r;
        pnode.a[1] = c;
        pnode.length = 1;
        pnode.score = 1;
        
        double pscore = 1;            
        for(j= 0; j< numchild; j++)
        {
            childid = flatparts[0].dyx[partids.a[i]].a[3+5*j];
            // define the children node
            struct Node cnode;
            // check its types
            int* omatched = find(o2list, numdyx, childid, 'e');;
            nummatched = sum(omatched, numdyx);
            int *bestk = (int*)malloc(numdyx*sizeof(int));
            //check its relative positions
            maxprob = 0;// for best match from a bunch of instances
            if(nummatched>0) 
            {
                 //then match relative positions
                 mean[0] = flatparts[0].dyx[partids.a[i]].a[4+5*j];
                 mean[1] = flatparts[0].dyx[partids.a[i]].a[5+5*j];
                 var[0] = flatparts[0].dyx[partids.a[i]].a[6+5*j];
                 var[1] = flatparts[0].dyx[partids.a[i]].a[7+5*j];
                 //find a best one in the current neighs
                 for(k=0; k<numdyx; k++)
                 {
                     if(omatched[k]==1)
                     {
                         double prob = (0.5*(dyx[2*k]- mean[0])*(dyx[2*k]- mean[0]))/(var[0]+4.5*pow(epsilon, numlayer-1)) + 0.5*((dyx[2*k+1]- mean[1])*(dyx[2*k+1]- mean[1]))/(var[1]+4.5*pow(epsilon,numlayer-1)) ;
                         //prob = exp(-prob);
                         if(maxprob<exp(-prob))
                         {
                             maxprob = exp(-prob);
                             bestk[j] = k;
                         }
                     }
                 }

                 if(maxprob > 0.5*pow(discount,numlayer-1))
                 { //it means a match
                    //int cy, cx, co2, radius;
                    pscore = pscore*maxprob;
                    cy = r + dyx[2*bestk[j]];
                    cx = c + dyx[2*bestk[j]+1];
                    co2 = o2list[bestk[j]];             
                    //get its children from the previous layer
                    if (numlayer>2)
                    {    
                        // find a neighborhood
                        radius = model[numlayer-2].radius;
                        reidx = searchneighs(cy, cx, pinstances, ninst, col, radius, numlayer-1);
                        pre_radius = 2;
                        if(numlayer>3)
                            pre_radius = model[numlayer-3].radius;
                        idx1 = searchneighs2(cy, cx, pinstances, ninst, col, pre_radius, numlayer-2, radjust);
                        // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
                        reidx = setdiff(reidx,idx1);
                        
                        int *y = (int*)malloc(sizeof(int)*reidx.length); 
                        // int x[idx.length];
                        int *x = (int*)malloc(sizeof(int)*reidx.length); 
                        int *iy = (int*)malloc(sizeof(int)*reidx.length); 
                        int *ix = (int*)malloc(sizeof(int)*reidx.length); 
                        int *co2list = (int*)malloc(sizeof(int)*reidx.length); 
                        int *cdyx = (int*)malloc(2*reidx.length*sizeof(int));
                        
                        for(int ii=0;ii< reidx.length; ii++)
                        {

                            y[ii] = pinstances[px2(reidx.a[ii], 1, ninst, col)];//pinstances[idx.a[i]][2];
                            x[ii] = pinstances[px2(reidx.a[ii], 2, ninst, col)];//pinstances[idx.a[i]][3]; 
                            co2list[ii] = pinstances[px2(reidx.a[ii], 0, ninst, col)];//pinstances[idx.a[i]][1];
                            iy[ii] = y[ii] - cy + radius;
                            ix[ii] = x[ii] -cx + radius;
                            cdyx[2*ii]= iy[ii]- radius;
                            cdyx[2*ii+1] = ix[ii]-radius; 

                        }
                        // matching it
                        cfound = topdownmatch(pinstances, ninst, col, co2, cy,cx, cdyx, co2list,reidx.length, 
                                model, numlayer-1, numOrient, numsub, radjust, discount, &cnode);
                    
                        free(y); y =NULL;
                        free(x); x = NULL;
                        free(iy); iy =NULL;
                        free(ix); ix = NULL;
                        free(co2list); co2list = NULL;
                        free(cdyx); cdyx = NULL;
                        //add here to speed up
                        if(!cfound)
                        {
                            if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                            if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                            
                            if(flag!=NULL)
                            {free(flag);flag=NULL;}
                            
                            //free some allocation/
                            /*if(indicator!=NULL)
                             {free(indicator);indicator = NULL;}
                             if(tscores!=NULL)
                             {free(tscores);tscores = NULL;}
                             if(partnodes!=NULL)
                             {free(partnodes);partnodes = NULL;}
                             if(tdlinks!=NULL)
                             {   
                                 free(tdlinks);
                                 tdlinks = NULL;
                             }
                             if(flatparts!=NULL)
                             {free(flatparts);flatparts = NULL;}
                             if(childlens!=NULL)
                             {free(childlens);childlens = NULL;}
                             return cfound;*/
                             break;
                        }
                    }
                    else{
                        cfound = true;
                        cnode.a[0] = cy;
                        cnode.a[1] = cx;
                        cnode.length = 1;
                        cnode.score = pscore;
                    }
                    
                    if(cfound){
                        flag[j] = 1;
                        pnode = cat(pnode, cnode);
                        pnode.score = pnode.score*cnode.score;
                    }
                 
                 }
                 //add here to speed up
                 else{
                        if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                        if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                        if(flag!=NULL)
                            {free(flag);flag=NULL;}
                        break;
                 }

            }
            if(omatched!=NULL)
            {free(omatched);omatched = NULL;}
            if(bestk!=NULL)
            {free(bestk);bestk = NULL;}


        }
        
        if(flag!=NULL)
        {
            //check the matching from its children
            if(sum(flag, numchild)==numchild){
                found =1;
                indicator[tindex] =1;
                tscores[tindex] = pnode.score;
                partnodes[tindex] = pnode;
            }
            free(flag);flag=NULL;
        }
    }
    
    /*how many two-subcomposition*/
    nummatched = 0;
    maxprob = 0;//remember the max
    int id=0; //remember id
    /* consider suppression here */
    for(i=0;i<partids.length; i++)
    {
        if(childlens[i]==2 && tscores[i]>0) //consider two children first
        {    
            nummatched = nummatched + 1;
            if(maxprob< tscores[i]) 
            {
                maxprob = tscores[i];
                id = i;
            }
        }
    }

    // only consider the 1-sub maximum score here
    if(nummatched==0)
        id= vmax(tscores, partids.length, &maxprob);
    /*copy result to output*/
    point->length =  partnodes[id].length;
    for(int jj=0; jj<2*point->length; jj++)
        point->a[jj] = partnodes[id].a[jj];
    point->compid = partids.a[id] +1;//compositional model id, +1 because it starts from 0
    point->score = maxprob;
    /*free some allocation*/
    if(indicator!=NULL)
    {free(indicator);indicator = NULL;}
    if(tscores!=NULL)
    {free(tscores);tscores = NULL;}
    if(partnodes!=NULL)
    {free(partnodes);partnodes = NULL;}
    if(tdlinks!=NULL)
    {   
        free(tdlinks);
        tdlinks = NULL;
    }
    if(flatparts!=NULL)
    {free(flatparts);flatparts = NULL;}
    if(childlens!=NULL)
    {free(childlens);childlens = NULL;}
    
    return found;
}


/*TOP DOWN MATCHING WITH ROTATION*/
bool topdownmatch_old(double *pinstances, int ninst, int col, int o1, int r, int c, int *dyx, int *o2list,int numdyx,
        Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point)
{   
    // initialize the value for return
    bool found, cfound;
    found= false; cfound = false;   

    struct array1d *tdlinks = (struct array1d*)malloc(numOrient*sizeof(struct array1d));
    struct pTypes *flatparts = mappingParts(model, numOrient, numlayer, numsub, rotation, tdlinks);

    //get the matching parts id
    struct array1d partids = tdlinks[o1-1];
    if(partids.length<1){
        found=false;
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        return found;
    } 
    int pDim =4;
    // define variables
    int i,j,k, pid;
    double epsilon = 1.5;
    double mean[2]; double var[2];
    int childid, numchild, anchorid;
    int cy, cx, co2, radius, pre_radius;//for its children
    struct array1d reidx, idx1;
    // initialize it
    int *indicator = (int*)malloc(partids.length*sizeof(int));
    double *tscores = (double*)malloc(partids.length*sizeof(double));
    int *childlens = (int*)malloc(partids.length*sizeof(int));
    for(i=0;i< partids.length; i++)
    {
        indicator[i] = 0;
        tscores[i] = 0;
        childlens[i] = 1;
    }
    //define the part nodes by the number of pids
    struct Node *partnodes = (struct Node*)malloc(partids.length*sizeof(struct Node));
    // check its neighborhood
    int tindex =-1;
    int nummatched = 0;
    double maxprob=0;
    #pragma omp parallel for 
    for(i=0; i<partids.length; i++)
    { // processing each part type
        tindex = tindex + 1;
        partnodes[i].length = 0; //initialize it
        pid = flatparts[0].dyx[partids.a[i]].a[0];//start from 1
        numchild = flatparts[0].dyx[partids.a[i]].a[2];
        childlens[i] = numchild;
        //indicator matching for child
        int* flag = (int*)malloc(numchild*sizeof(int));
        initialize(flag, numchild);
        //define the current node from this input point
        struct Node pnode;
        pnode.a[0] = r;
        pnode.a[1] = c;        
        pnode.a[2] = pid;//pid
        pnode.a[3] = numlayer;
        pnode.length = 1;
        pnode.score = 1;
        
        double pscore = 1;            
        for(j= 0; j< numchild; j++)
        {
            //model[numlayer-1].ptypes_[numsub-1].dyx[pid].a[0];
            childid = flatparts[0].dyx[partids.a[i]].a[3+5*j];
            // define the children node
            struct Node cnode;
            // check its types
            int* omatched = find(o2list, numdyx, childid, 'e');;
            nummatched = sum(omatched, numdyx);
            int *bestk = (int*)malloc(numchild*sizeof(int));
            //check its relative positions
            maxprob = 0;// for best match from a bunch of instances
            if(nummatched>0) 
            {
                 //then match relative positions
                 mean[0] = flatparts[0].dyx[partids.a[i]].a[4+5*j];
                 mean[1] = flatparts[0].dyx[partids.a[i]].a[5+5*j];
                 var[0] = flatparts[0].dyx[partids.a[i]].a[6+5*j];
                 var[1] = flatparts[0].dyx[partids.a[i]].a[7+5*j];
                 //find a best one in the current neighs
                 for(k=0; k<numdyx; k++)
                 {
                     if(omatched[k]==1)
                     {
                         double prob = (0.5*(dyx[2*k]- mean[0])*(dyx[2*k]- mean[0]))/(var[0]+4.5*pow(epsilon, numlayer-1)) + 0.5*((dyx[2*k+1]- mean[1])*(dyx[2*k+1]- mean[1]))/(var[1]+4.5*pow(epsilon,numlayer-1)) ;
                         //prob = exp(-prob);
                         if(maxprob<exp(-prob))
                         {
                             maxprob = exp(-prob);
                             bestk[j] = k;
                         }
                     }
                 }

                 if(maxprob > 0.5*pow(discount,numlayer-1))
                 { //it means a match
                    //int cy, cx, co2, radius;
                    pscore = pscore*maxprob;
                    cy = r + dyx[2*bestk[j]];
                    cx = c + dyx[2*bestk[j]+1];
                    co2 = o2list[bestk[j]];             
                    //get its children from the previous layer
                    if (numlayer>2)
                    {    
                        // find a neighborhood
                        radius = model[numlayer-2].radius;
                        reidx = searchneighs(cy, cx, pinstances, ninst, col, radius, numlayer-1);
                        pre_radius = 2;
                        if(numlayer>3)
                            pre_radius = model[numlayer-3].radius;
                        idx1 = searchneighs2(cy, cx, pinstances, ninst, col, pre_radius, numlayer-2, radjust);
                        // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
                        reidx = setdiff(reidx,idx1);
                        
                        int *y = (int*)malloc(sizeof(int)*reidx.length); 
                        // int x[idx.length];
                        int *x = (int*)malloc(sizeof(int)*reidx.length); 
                        int *iy = (int*)malloc(sizeof(int)*reidx.length); 
                        int *ix = (int*)malloc(sizeof(int)*reidx.length); 
                        int *co2list = (int*)malloc(sizeof(int)*reidx.length); 
                        int *cdyx = (int*)malloc(2*reidx.length*sizeof(int));
                        
                        for(int ii=0;ii< reidx.length; ii++)
                        {

                            y[ii] = pinstances[px2(reidx.a[ii], 1, ninst, col)];//pinstances[idx.a[i]][2];
                            x[ii] = pinstances[px2(reidx.a[ii], 2, ninst, col)];//pinstances[idx.a[i]][3]; 
                            co2list[ii] = pinstances[px2(reidx.a[ii], 0, ninst, col)];//pinstances[idx.a[i]][1];
                            iy[ii] = y[ii] - cy + radius;
                            ix[ii] = x[ii] -cx + radius;
                            cdyx[2*ii]= iy[ii]- radius;
                            cdyx[2*ii+1] = ix[ii]-radius; 

                        }
                        // matching it
                        cfound = topdownmatch(pinstances, ninst, col, co2, cy,cx, cdyx, co2list,reidx.length, 
                                model, numlayer-1, numOrient, numsub, radjust, discount, rotation, &cnode);
                    
                        free(y); y =NULL;
                        free(x); x = NULL;
                        free(iy); iy =NULL;
                        free(ix); ix = NULL;
                        free(co2list); co2list = NULL;
                        free(cdyx); cdyx = NULL;
                        //add here to speed up
                        if(!cfound)
                        {
                            if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                            if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                            
                            if(flag!=NULL)
                            {free(flag);flag=NULL;}
                            
                            //free some allocation/
                            /*if(indicator!=NULL)
                             {free(indicator);indicator = NULL;}
                             if(tscores!=NULL)
                             {free(tscores);tscores = NULL;}
                             if(partnodes!=NULL)
                             {free(partnodes);partnodes = NULL;}
                             if(tdlinks!=NULL)
                             {   
                                 free(tdlinks);
                                 tdlinks = NULL;
                             }
                             if(flatparts!=NULL)
                             {free(flatparts);flatparts = NULL;}
                             if(childlens!=NULL)
                             {free(childlens);childlens = NULL;}
                             return cfound;*/
                             break;
                        }/*
                        else
                            {
                                // add current matching points to cnode 
                                cnode.a[cnode.length*pDim+0] = cy;
                                cnode.a[cnode.length*pDim+1] = cx;
                                cnode.a[cnode.length*pDim+2]= model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                                cnode.a[cnode.length*pDim+3] = numlayer-1;//its child from previous layer
                                cnode.length = cnode.length + 1;
                                
                            }*/
                    }
                    else{
                        
                        cfound = true;
                        cnode.a[0] = cy;
                        cnode.a[1] = cx;
                        cnode.a[2]= model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                        cnode.a[3] = numlayer-1;//its child from previous layer
                        cnode.length = 1;
                        cnode.score = maxprob;
                    }
                    
                    if(cfound){
                        flag[j] = 1;
                        pnode = cat(pnode, cnode, pDim);
                        pnode.score = pnode.score*cnode.score;
                    }
                 
                 }
                 //add here to speed up
                 else{
                        if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                        if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                        if(flag!=NULL)
                            {free(flag);flag=NULL;}
                        break;
                 }

            }
            if(omatched!=NULL)
            {free(omatched);omatched = NULL;}
            if(bestk!=NULL)
            {free(bestk);bestk = NULL;}


        }
        
        if(flag!=NULL)
        {
            //check the matching from its children
            if(sum(flag, numchild)==numchild){
                found =1;
                indicator[tindex] =1;
                tscores[tindex] = pnode.score;
                partnodes[tindex] = pnode;
            }
            //else found = 0;
            free(flag);flag=NULL;
        }
        //else
        //    found = 0;
    }
    
    /*how many two-subcomposition*/
    nummatched = 0;
    maxprob = 0;//remember the max
    int id=0; //remember id
    /* consider suppression here */
    for(i=0;i<partids.length; i++)
    {
        if(childlens[i]==2 && tscores[i]>0) //consider two children first
        {    
            nummatched = nummatched + 1;
            if(maxprob< tscores[i]) 
            {
                maxprob = tscores[i];
                id = i;
            }
        }
    }

    // only consider the 1-sub maximum score here
    if(nummatched==0)
        id= vmax(tscores, partids.length, &maxprob);
    /*copy result to output*/
    point->length =  partnodes[id].length;
    for(int jj=0; jj<pDim*point->length; jj++)
        point->a[jj] = partnodes[id].a[jj];
    point->compid = partids.a[id] +1;//compositional model id, +1 because it starts from 0
    point->score = maxprob;
    /*free some allocation*/
    if(indicator!=NULL)
    {free(indicator);indicator = NULL;}
    if(tscores!=NULL)
    {free(tscores);tscores = NULL;}
    if(partnodes!=NULL)
    {free(partnodes);partnodes = NULL;}
    if(tdlinks!=NULL)
    {   
        free(tdlinks);
        tdlinks = NULL;
    }
    if(flatparts!=NULL)
    {free(flatparts);flatparts = NULL;}
    if(childlens!=NULL)
    {free(childlens);childlens = NULL;}
    
    return found;
}


/*TOP DOWN MATCHING WITH ROTATION*/
bool topdownmatch(double *pinstances, int ninst, int col, int o1, int r, int c, int *dyx, int *o2list,int numdyx,
        Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point)
{   
    // initialize the value for return
    bool found, cfound;
    found= false; cfound = false;   

    struct array1d *tdlinks = (struct array1d*)malloc(numOrient*sizeof(struct array1d));
    struct pTypes *flatparts = mappingParts(model, numOrient, numlayer, numsub, rotation, tdlinks);
    if(flatparts==NULL)
    {   
        
        point->length =  1;
        point->a[0]= r;
        point->a[1] =c;
        point->a[2] = o1;
        point->a[3] = numlayer;
        point->compid = o1;//compositional model id, +1 because it starts from 0
        point->score = 1;
        
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        return true;
    }
    //get the matching parts id
    struct array1d partids = tdlinks[o1-1];
    if(partids.length<1){
        found=false;
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        return found;
    } 
    int pDim =4;
    // define variables
    int i,j,k, pid;
    double epsilon = 1.5;
    double mean[2]; double var[2];
    int childid, numchild, anchorid;
    int cy, cx, co2, radius, pre_radius;//for its children
    struct array1d reidx, idx1;
    // initialize it
    int *indicator = (int*)malloc(partids.length*sizeof(int));
    double *tscores = (double*)malloc(partids.length*sizeof(double));
    int *childlens = (int*)malloc(partids.length*sizeof(int));
    for(i=0;i< partids.length; i++)
    {
        indicator[i] = 0;
        tscores[i] = 0;
        childlens[i] = 1;
    }
    //define the part nodes by the number of pids
    struct Node *partnodes = (struct Node*)malloc(partids.length*sizeof(struct Node));
    // check its neighborhood
    int tindex =-1;
    int nummatched = 0;
    double maxprob=0;
    #pragma omp parallel for 
    for(i=0; i<partids.length; i++)
    { // processing each part type
        tindex = tindex + 1;
        partnodes[i].length = 0; //initialize it
        pid = flatparts[0].dyx[partids.a[i]].a[0];//start from 1
        numchild = flatparts[0].dyx[partids.a[i]].a[2];
        childlens[i] = numchild;
        //indicator matching for child
        int* flag = (int*)malloc(numchild*sizeof(int));
        initialize(flag, numchild);
        //define the current node from this input point
        struct Node pnode;
        pnode.a[0] = r;
        pnode.a[1] = c;        
        pnode.a[2] = pid;//pid
        pnode.a[3] = numlayer;
        pnode.length = 1;
        pnode.score = 1;
        
        double pscore = 1;            
        for(j= 0; j< numchild; j++)
        {
            //model[numlayer-1].ptypes_[numsub-1].dyx[pid].a[0];
            childid = flatparts[0].dyx[partids.a[i]].a[3+5*j];
            // define the children node
            struct Node cnode;
            // check its types
            int* omatched = find(o2list, numdyx, childid, 'e');;
            nummatched = sum(omatched, numdyx);
            int *bestk = (int*)malloc(numchild*sizeof(int));
            //check its relative positions
            maxprob = 0;// for best match from a bunch of instances
            if(nummatched>0) 
            {
                 //then match relative positions
                 mean[0] = flatparts[0].dyx[partids.a[i]].a[4+5*j];
                 mean[1] = flatparts[0].dyx[partids.a[i]].a[5+5*j];
                 var[0] = flatparts[0].dyx[partids.a[i]].a[6+5*j];
                 var[1] = flatparts[0].dyx[partids.a[i]].a[7+5*j];
                 //find a best one in the current neighs
                 for(k=0; k<numdyx; k++)
                 {
                     if(omatched[k]==1)
                     {
                         double prob = (0.5*(dyx[2*k]- mean[0])*(dyx[2*k]- mean[0]))/(var[0]+4.5*pow(epsilon, numlayer-1)) + 0.5*((dyx[2*k+1]- mean[1])*(dyx[2*k+1]- mean[1]))/(var[1]+4.5*pow(epsilon,numlayer-1)) ;
                         //prob = exp(-prob);
                         if(maxprob<exp(-prob))
                         {
                             maxprob = exp(-prob);
                             bestk[j] = k;
                         }
                     }
                 }

                 if(maxprob > 0.5*pow(discount,numlayer-1))
                 { //it means a match
                    //int cy, cx, co2, radius;
                    pscore = pscore*maxprob;
                    cy = r + dyx[2*bestk[j]];
                    cx = c + dyx[2*bestk[j]+1];
                    co2 = o2list[bestk[j]];             
                    //get its children from the previous layer
                    if (numlayer>1)
                    {    
                        // find a neighborhood
                        radius = model[numlayer-2].radius;
                        reidx = searchneighs(cy, cx, pinstances, ninst, col, radius, numlayer-1);
                        pre_radius = 2;
                        if(numlayer>3)
                            pre_radius = model[numlayer-3].radius;
                        idx1 = searchneighs2(cy, cx, pinstances, ninst, col, pre_radius, numlayer-2, radjust);
                        // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
                        reidx = setdiff(reidx,idx1);
                        
                        int *y = (int*)malloc(sizeof(int)*reidx.length); 
                        // int x[idx.length];
                        int *x = (int*)malloc(sizeof(int)*reidx.length); 
                        int *iy = (int*)malloc(sizeof(int)*reidx.length); 
                        int *ix = (int*)malloc(sizeof(int)*reidx.length); 
                        int *co2list = (int*)malloc(sizeof(int)*reidx.length); 
                        int *cdyx = (int*)malloc(2*reidx.length*sizeof(int));
                        
                        for(int ii=0;ii< reidx.length; ii++)
                        {

                            y[ii] = pinstances[px2(reidx.a[ii], 1, ninst, col)];//pinstances[idx.a[i]][2];
                            x[ii] = pinstances[px2(reidx.a[ii], 2, ninst, col)];//pinstances[idx.a[i]][3]; 
                            co2list[ii] = pinstances[px2(reidx.a[ii], 0, ninst, col)];//pinstances[idx.a[i]][1];
                            iy[ii] = y[ii] - cy + radius;
                            ix[ii] = x[ii] -cx + radius;
                            cdyx[2*ii]= iy[ii]- radius;
                            cdyx[2*ii+1] = ix[ii]-radius; 

                        }
                        // matching it
                        cfound = topdownmatch(pinstances, ninst, col, co2, cy,cx, cdyx, co2list,reidx.length, 
                                model, numlayer-1, numOrient, numsub, radjust, discount, rotation, &cnode);
                    
                        free(y); y =NULL;
                        free(x); x = NULL;
                        free(iy); iy =NULL;
                        free(ix); ix = NULL;
                        free(co2list); co2list = NULL;
                        free(cdyx); cdyx = NULL;
                        //add here to speed up
                        if(!cfound)
                        {
                            if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                            if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                            
                            if(flag!=NULL)
                            {free(flag);flag=NULL;}
                            
                            //free some allocation/
                            /*if(indicator!=NULL)
                             {free(indicator);indicator = NULL;}
                             if(tscores!=NULL)
                             {free(tscores);tscores = NULL;}
                             if(partnodes!=NULL)
                             {free(partnodes);partnodes = NULL;}
                             if(tdlinks!=NULL)
                             {   
                                 free(tdlinks);
                                 tdlinks = NULL;
                             }
                             if(flatparts!=NULL)
                             {free(flatparts);flatparts = NULL;}
                             if(childlens!=NULL)
                             {free(childlens);childlens = NULL;}
                             return cfound;*/
                             break;
                        }
                        else
                            {
                                // add current matching points to cnode 
                                cnode.a[cnode.length*pDim+0] = cy;
                                cnode.a[cnode.length*pDim+1] = cx;
                                cnode.a[cnode.length*pDim+2]= model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                                cnode.a[cnode.length*pDim+3] = numlayer-1;//its child from previous layer
                                cnode.length = cnode.length + 1;
                                
                            }
                    }
                    else{
                        /*
                        cfound = true;
                        cnode.a[0] = cy;
                        cnode.a[1] = cx;
                        cnode.a[2]= model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                        cnode.a[3] = numlayer-1;//its child from previous layer
                        cnode.length = 1;
                        cnode.score = maxprob;
                        */
                    }
                    
                    if(cfound){
                        flag[j] = 1;
                        pnode = cat(pnode, cnode, pDim);
                        pnode.score = pnode.score*cnode.score;
                    }
                 
                 }
                 //add here to speed up
                 else{
                        if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                        if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                        if(flag!=NULL)
                            {free(flag);flag=NULL;}
                        break;
                 }

            }
            if(omatched!=NULL)
            {free(omatched);omatched = NULL;}
            if(bestk!=NULL)
            {free(bestk);bestk = NULL;}


        }
        
        if(flag!=NULL)
        {
            //check the matching from its children
            if(sum(flag, numchild)==numchild){
                found =1;
                indicator[tindex] =1;
                tscores[tindex] = pnode.score;
                partnodes[tindex] = pnode;
            }
            //else found = 0;
            free(flag);flag=NULL;
        }
        //else
        //    found = 0;
    }
    
    /*how many two-subcomposition*/
    nummatched = 0;
    maxprob = 0;//remember the max
    int id=0; //remember id
    /* consider suppression here */
    for(i=0;i<partids.length; i++)
    {
        if(childlens[i]==2 && tscores[i]>0) //consider two children first
        {    
            nummatched = nummatched + 1;
            if(maxprob< tscores[i]) 
            {
                maxprob = tscores[i];
                id = i;
            }
        }
    }

    // only consider the 1-sub maximum score here
    if(nummatched==0)
        id= vmax(tscores, partids.length, &maxprob);
    /*copy result to output*/
    point->length =  partnodes[id].length;
    for(int jj=0; jj<pDim*point->length; jj++)
        point->a[jj] = partnodes[id].a[jj];
    point->compid = partids.a[id] +1;//compositional model id, +1 because it starts from 0
    point->score = maxprob;
    /*free some allocation*/
    if(indicator!=NULL)
    {free(indicator);indicator = NULL;}
    if(tscores!=NULL)
    {free(tscores);tscores = NULL;}
    if(partnodes!=NULL)
    {free(partnodes);partnodes = NULL;}
    if(tdlinks!=NULL)
    {   
        free(tdlinks);
        tdlinks = NULL;
    }
    if(flatparts!=NULL)
    {free(flatparts);flatparts = NULL;}
    if(childlens!=NULL)
    {free(childlens);childlens = NULL;}
    
    return found;
}

/*TOP DOWN MATCHING WITH ROTATION, deleting matched point in the next recursive step*/
bool topdownmatch_fast(double *pinstances, int ninst, int col, int *pdist, int iinst, struct array1d idx, bool *visited,
        Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point)
{   
    // initialize the value for return
    bool found, cfound;
    found= false; 
    cfound = false; 
    
    int o1, r,c;
    int *y,*x,*iy,*ix, *dyx,*o2list;
    //initialize the location and orientation
    r = pinstances[px2(iinst, 1, ninst, pDim)];
    c = pinstances[px2(iinst, 2, ninst, pDim)];
    o1 = pinstances[px2(iinst, 0, ninst, pDim)];
    //delete selected points by setting flag
    pinstances[px2(iinst, 3, ninst, pDim)] =0;
    //printf("r=%d, c=%d, iinst=%d ", r, c, iinst);
    struct array1d *tdlinks = (struct array1d*)malloc(numOrient*sizeof(struct array1d));
    struct pTypes *flatparts = mappingParts(model, numOrient, numlayer, numsub, rotation, tdlinks);
    /*
    if(flatparts==NULL)
    {   
        
        point->length =  1;
        point->a[0]= r;
        point->a[1] =c;
        point->a[2] = o1;
        point->a[3] = numlayer;
        point->compid = o1;//compositional model id
        point->score = 1;
        
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        visited[iinst] = true;
        return true;
    }*/
    //get the matching parts id
    struct array1d partids = tdlinks[o1-1];
    if(partids.length<1){
        found=false;
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        return found;
    } 
    int pDim =4;
    // define variables
    int i,j,k, pid;
    double epsilon = 1.5;
    double mean[2]; double var[2];
    int childid, numchild, anchorid;
    int cy, cx, co2, radius, pre_radius;//for its children
    struct array1d reidx;
    
    /*
    // show links
    printf("---------------- begin -----------\n");
    for(i=0;i<numOrient;i++)
    {
        printf("the %d-th orient with length %d is: \n", i,tdlinks[i].length);
        for(j=0;j<tdlinks[i].length;j++)
            printf(" %d", tdlinks[i].a[j]);
        printf("end\n", i);
    }
    printf("---------------- end --------------\n");
    */
    
    //get its neighbors
    if(idx.length>0)
    {
        y = (int*)malloc(sizeof(int)*idx.length); 
        //int x[idx.length];
        x = (int*)malloc(sizeof(int)*idx.length); 
        iy = (int*)malloc(sizeof(int)*idx.length); 
        ix = (int*)malloc(sizeof(int)*idx.length); 
        o2list = (int*)malloc(sizeof(int)*idx.length); 
        dyx = (int*)malloc(2*idx.length*sizeof(int));
        for(i=0;i< idx.length; i++)
        {

            y[i] = (int)pinstances[px2(idx.a[i], 1, ninst, pDim)];//pinstances[idx.a[i]][2];
            x[i] = (int)pinstances[px2(idx.a[i], 2, ninst, pDim)];//pinstances[idx.a[i]][3]; 
            o2list[i] = (int)pinstances[px2(idx.a[i], 0, ninst, pDim)];//pinstances[idx.a[i]][1];
            iy[i] = y[i] - r + radius;
            ix[i] = x[i] -c + radius;
            dyx[2*i]= iy[i]- radius;
            dyx[2*i+1] = ix[i]-radius; 

        }

    }
    else{
        
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        
        return found;
    }
    // initialize the structure and score for each matched part id
    int *indicator = (int*)malloc(partids.length*sizeof(int));
    double *tscores = (double*)malloc(partids.length*sizeof(double));
    int *childlens = (int*)malloc(partids.length*sizeof(int));
    for(i=0;i< partids.length; i++)
    {
        indicator[i] = 0;
        tscores[i] = 0;
        childlens[i] = 1;
        //flatparts[0].dyx[partids.a[i]].a[2];
    }
    //define the part nodes by the number of pids
    struct Node *partnodes = (struct Node*)malloc(partids.length*sizeof(struct Node));
    // processing each part type
    int tindex =-1;
    int nummatched = 0;
    double maxprob=0;
    #pragma omp parallel for 
    for(i=0; i<partids.length; i++)
    { 
        tindex = tindex + 1;
        partnodes[i].length = 0; //initialize it
        pid = flatparts[0].dyx[partids.a[i]].a[0];//start from 1
        numchild = flatparts[0].dyx[partids.a[i]].a[2];
        childlens[i] = numchild;
        //printf("numchild=%d, numlayer=%d\n", numchild, numlayer);
        //define the current node from this input point
        
        struct Node pnode;
        pnode.a[0] = r;
        pnode.a[1] = c;        
        pnode.a[2] = pid;//partids.a[i];
        pnode.a[3] = numlayer;
        pnode.length = 1;
        pnode.score = 1;
        
        double pscore = 1;   
        
        //indicator matching for child
        int* flag = (int*)malloc(numchild*sizeof(int));
        initialize(flag, numchild);
        int *bestk = (int*)malloc(numchild*sizeof(int));   
        
        //mask current point
        visited[iinst] = true;
        // check subparts with its neighborhood
        for(j= 0; j< numchild; j++)
        {
            //model[numlayer-1].ptypes_[numsub-1].dyx[pid].a[0];
            childid = flatparts[0].dyx[partids.a[i]].a[3+5*j];
            // define the children node
            struct Node cnode;
            cnode.length = 0;
            // check its types
            int* omatched = find(o2list, idx.length, childid, 'e');;
            nummatched = sum(omatched, idx.length);
            //check its relative positions
            maxprob = 0;// for best match from a bunch of instances
            if(nummatched>0) 
            {
                 //then match relative positions
                 mean[0] = flatparts[0].dyx[partids.a[i]].a[4+5*j];
                 mean[1] = flatparts[0].dyx[partids.a[i]].a[5+5*j];
                 var[0] = flatparts[0].dyx[partids.a[i]].a[6+5*j];
                 var[1] = flatparts[0].dyx[partids.a[i]].a[7+5*j];
                 //find a best one in the current neighs
                 for(k=0; k<idx.length; k++)
                 {
                     if(omatched[k]==1)
                     {
                         double prob = (0.5*(dyx[2*k]- mean[0])*(dyx[2*k]- mean[0]))/(var[0]+4.5*pow(epsilon, numlayer-1)) + 0.5*((dyx[2*k+1]- mean[1])*(dyx[2*k+1]- mean[1]))/(var[1]+4.5*pow(epsilon,numlayer-1)) ;
                         //prob = exp(-prob);
                         if(maxprob<exp(-prob))
                         {
                             maxprob = exp(-prob);
                             bestk[j] = k;
                         }
                     }
                 }

                 if(maxprob > 0.5*pow(discount,numlayer-1))
                 { //it means a match

                    //get the score and its subpart locations int cy, cx, co2, radius;
                    pscore = pscore*maxprob;
                    cy = r + dyx[2*bestk[j]];
                    cx = c + dyx[2*bestk[j]+1];
                    co2 = o2list[bestk[j]];             
                    //get its children from the previous layer
                    if (numlayer>2)
                    {    
                        // find a neighborhood
                        radius = model[numlayer-2].radius; 
                        pre_radius = 2;
                        if(numlayer>3)
                            pre_radius = model[numlayer-3].radius;// - radjust[numlayer-2];
                        reidx = searchneighs(idx.a[bestk[j]], pdist, ninst, radius, pre_radius, visited);
                        //printf("current id: %d, the number of neighbors: %d\n", idx.a[bestk[j]], reidx.length);
                        if(reidx.length>0 && !visited[idx.a[bestk[j]]])
                            // matching it
                            cfound = topdownmatch_fast(pinstances, ninst, pDim, pdist, idx.a[bestk[j]], reidx, visited, model, numlayer-1, numOrient, 
                        numsub, radjust, discount,rotation, &cnode);
                        
                        if(!cfound)
                        {
                            if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                            if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                            if(flag!=NULL)
                            {free(flag);flag=NULL;}
                            
                             break;
                        }/*
                        else
                            {
                                
                                visited[idx.a[bestk[j]]] = true;
                                // add current matching points to cnode 
                                cnode.a[cnode.length*pDim+0] = cy;
                                cnode.a[cnode.length*pDim+1] = cx;
                                cnode.a[cnode.length*pDim+2]= model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                                cnode.a[cnode.length*pDim+3] = numlayer-1;//its child from previous layer
                                cnode.length = cnode.length + 1;
                                
                            }
                          */
                    }
                    else{
                        
                        cfound = true;
                        
                        //visited[idx.a[bestk[j]]] = true;
                        cnode.a[0] = cy;
                        cnode.a[1] = cx;
                        cnode.a[2]= model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                        cnode.a[3] = numlayer-1;//its child from previous layer
                        cnode.length = 1;
                        cnode.score = maxprob; 
                        // other considerations
                        // cnode.length = 0;
                        /*
                        cnode.a[cnode.length*pDim+0] = cy;
                        cnode.a[cnode.length*pDim+1] = cx;
                        cnode.a[cnode.length*pDim+2]= model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                        cnode.a[cnode.length*pDim+3] = numlayer-1;//its child from previous layer
                        cnode.length = cnode.length + 1;
                        */
                        
                        
                    }
                    
                    if(cfound){
                        flag[j] = 1;
                        /*                        
                        // add children here
                        pnode.a[pnode.length*pDim+0] = cy;
                        pnode.a[pnode.length*pDim+1] = cx;        
                        pnode.a[pnode.length*pDim+2] = model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];//pid
                        pnode.a[pnode.length*pDim+3] = numlayer-1;
                        pnode.length = pnode.length+1;
                        
                        for(int ci=0; ci<cnode.length;ci++)
                        {
                            printf("cy=%d, cx= %d, ", cnode.a[ci*pDim], cnode.a[ci*pDim+1]);
                        }
                        printf("\n");*/
                        pnode = cat(pnode, cnode, pDim);
                        // printf("the length with one children: %d\n", pnode.length);
                        //cat(&pnode, &cnode, pDim);
                        pnode.score = pnode.score*cnode.score;
                    }
                 
                 }
                 //add here to speed up
                 else{
                        if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                        flag[j] =0;
                        break;
                 }

            }
            if(omatched!=NULL)
            {free(omatched);omatched = NULL;}

        }
        
        if(flag!=NULL)
        {
            //check the matching from its children
            if(sum(flag, numchild)==numchild){
                found =1;
                indicator[tindex] =1;
                tscores[tindex] = pnode.score;
                partnodes[tindex] = pnode;
                // mask it
                visited[iinst] = true;
                // mask its children
                for(j =0; j<numchild; j++)
                    visited[idx.a[bestk[j]]] = true;
            }
            //else found = 0;
            free(flag);flag=NULL;
        }
        //else
        //    found = 0;
        if(bestk!=NULL)
        {free(bestk);bestk = NULL;} 
        
    }
    
    //free neighbors
    if(y!=NULL)
    {   free(y); y=NULL;}
    if(x!=NULL)
    {   free(x); x= NULL;}
    if(iy!=NULL)
    {   free(iy); iy= NULL;}
    if(ix!= NULL)
    {   free(ix); ix = NULL;}
    if(o2list!= NULL)
    {   free(o2list); o2list = NULL;}
    if(dyx!= NULL)
    {   free(dyx);dyx = NULL;}
    
    
    /*check how many two-subcomposition matched*/
    nummatched = 0;
    maxprob = 0;//remember the max
    int id=0; //remember id
    /* consider suppression here */
    for(i=0;i<partids.length; i++)
    {
        if(childlens[i]==2 && tscores[i]>0) //consider two children first
        {    
            nummatched = nummatched + 1;
            if(maxprob< tscores[i]) 
            {
                maxprob = tscores[i];
                id = i;
            }
        }
    }

    // only consider the 1-sub maximum score here
    if(nummatched==0)
        id= vmax(tscores, partids.length, &maxprob);
    /*copy result to output*/
    point->length =  partnodes[id].length;
    for(int jj=0; jj<pDim*point->length; jj++)
        point->a[jj] = partnodes[id].a[jj];
    point->compid = partids.a[id] +1;//compositional model id, +1 because it starts from 0
    point->score = maxprob;
    /*free some allocation*/
    if(indicator!=NULL)
    {free(indicator);indicator = NULL;}
    if(tscores!=NULL)
    {free(tscores);tscores = NULL;}
    if(partnodes!=NULL)
    {free(partnodes);partnodes = NULL;}
    if(tdlinks!=NULL)
    {   
        free(tdlinks);
        tdlinks = NULL;
    }
    if(flatparts!=NULL)
    {free(flatparts);flatparts = NULL;}
    if(childlens!=NULL)
    {free(childlens);childlens = NULL;}
    
    return found;
}


/*TOP DOWN MATCHING WITH ROTATION, deleting matched point in the next recursive step*/
bool topdownmatch_full(double *pinstances, int ninst, int col, int *pdist, int iinst, struct array1d idx, bool *visited,
        Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point)
{   
    // initialize the value for return
    bool found, afound, cfound;
    found= false; cfound = false; 
    afound = false;
    int o1, r,c;
    int *y,*x,*iy,*ix, *dyx,*o2list;
    //initialize the location and orientation
    r = pinstances[px2(iinst, 1, ninst, pDim)];
    c = pinstances[px2(iinst, 2, ninst, pDim)];
    o1 = pinstances[px2(iinst, 0, ninst, pDim)];
    //delete selected points by setting flag
    pinstances[px2(iinst, 3, ninst, pDim)] =0;
    
    struct array1d *tdlinks = (struct array1d*)malloc(numOrient*sizeof(struct array1d));
    struct pTypes *flatparts = mappingParts(model, numOrient, numlayer, numsub, rotation, tdlinks);

    //get the matching parts id
    struct array1d partids = tdlinks[o1-1];
    if(partids.length<1){
        found=false;
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        return found;
    } 
    int pDim =4;
    // define variables
    int i,j,k;
    double epsilon = 1.5;
    double mean[2]; double var[2];
    int pid, childid, numchild, anchorid;
    int cy, cx, co2, radius, pre_radius;//for its children
    struct array1d reidx;
    
    radius = model[numlayer-1].radius; 
    //get its neighbors
    if(idx.length>0)
    {
        y = (int*)malloc(sizeof(int)*idx.length); 
        //int x[idx.length];
        x = (int*)malloc(sizeof(int)*idx.length); 
        iy = (int*)malloc(sizeof(int)*idx.length); 
        ix = (int*)malloc(sizeof(int)*idx.length); 
        o2list = (int*)malloc(sizeof(int)*idx.length); 
        dyx = (int*)malloc(2*idx.length*sizeof(int));
        for(i=0;i< idx.length; i++)
        {

            y[i] = (int)pinstances[px2(idx.a[i], 1, ninst, pDim)];//pinstances[idx.a[i]][2];
            x[i] = (int)pinstances[px2(idx.a[i], 2, ninst, pDim)];//pinstances[idx.a[i]][3]; 
            o2list[i] = (int)pinstances[px2(idx.a[i], 0, ninst, pDim)];//pinstances[idx.a[i]][1];
            iy[i] = y[i] - r + radius;
            ix[i] = x[i] -c + radius;
            dyx[2*i]= iy[i]- radius;
            dyx[2*i+1] = ix[i]-radius; 

        }

    }
    else{
        
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        
        return found;
    }
    // initialize the structure and score for each matched part id
    int *indicator = (int*)malloc(partids.length*sizeof(int));
    double *tscores = (double*)malloc(partids.length*sizeof(double));
    int *childlens = (int*)malloc(partids.length*sizeof(int));
    for(i=0;i< partids.length; i++)
    {
        indicator[i] = 0;
        tscores[i] = 0;
        childlens[i] = 1;
    }
    //define the part nodes by the number of pids
    struct Node *partnodes = (struct Node*)malloc(partids.length*sizeof(struct Node));
    // processing each part type
    int tindex =-1;
    int nummatched = 0;
    double maxprob=0;
    #pragma omp parallel for 
    for(i=0; i<partids.length; i++)
    { 
        tindex = tindex + 1;
        partnodes[i].length = 0; //initialize it
        pid = flatparts[0].dyx[partids.a[i]].a[0];//start from 1
        numchild = flatparts[0].dyx[partids.a[i]].a[2];
        childlens[i] = numchild;
        //printf("numchild=%d, numlayer=%d\n", numchild, numlayer);
        //define the current node from this input point
        struct Node pnode;
        /*pnode.a[0] = r;
        pnode.a[1] = c;        
        pnode.a[2] = pid;//partids.a[i];
        pnode.a[3] = numlayer;
        pnode.length = 1;
        pnode.score = 1;
        */
        pnode.length = 0;
        pnode.score = 1;
        double pscore = 1;   
        
        
        /*------------------------------------------------------------
         * matching for anchor point
         *------------------------------------------------------------
         */
        
        struct Node anode;
        anode.length = 0; //initialize it
        anode.score = 1;
        anchorid = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[1];
        if (numlayer>2)
        {    
            radius = model[numlayer-2].radius; 
            pre_radius = 2;
            if(numlayer>2)
                pre_radius = model[numlayer-3].radius;// - radjust[numlayer-2];
            reidx = searchneighs(iinst, pdist, ninst, radius, pre_radius, visited);
            // printf("current id: %d, the number of neighbors: %d", idx.a[bestk[j]], reidx.length);
            if(reidx.length>0 && !visited[iinst])
                afound = topdownmatch_full(pinstances, ninst, pDim, pdist, iinst, reidx, visited, model, numlayer-1, numOrient, 
                        numsub, radjust, discount,rotation, &anode);
        }
        else
        {

            afound = true;
            
            //visited[idx.a[bestk[j]]] = true;
            anode.a[0] = r;
            anode.a[1] = c;
            anode.a[2]= anchorid;//model[numlayer-1].ptypes_[numsub-1].dyx[pid].a[3+5*j];
            anode.a[3] = numlayer-1;//its child from previous layer
            anode.length = 1;
            anode.score = pscore;            
        }
        if(afound)
        {    
             //cat(&pnode, &anode, pDim);
             //pnode->score = pnode->score*cnode.score;
             pnode = cat(pnode, anode, pDim);
             pnode.score = pnode.score*anode.score;
        }
        
        
         /*------------------------------------------------------------
         * matching for childrens
         *------------------------------------------------------------
         */
        
        //indicator matching for child
        int* flag = (int*)malloc(numchild*sizeof(int));
        int *bestk = (int*)malloc(numchild*sizeof(int));   
        //initialize it
        initialize(flag, numchild);
        //mask current point
        visited[iinst] = true;
        // check subparts with its neighborhood
        for(j= 0; j< numchild; j++)
        {
            //model[numlayer-1].ptypes_[numsub-1].dyx[pid].a[0];
            childid = flatparts[0].dyx[partids.a[i]].a[3+5*j];
            // define the children node
            struct Node cnode;
            cnode.length = 0;
            // check its types
            int* omatched = find(o2list, idx.length, childid, 'e');;
            nummatched = sum(omatched, idx.length);
            //check its relative positions
            maxprob = 0;// for best match from a bunch of instances
            if(nummatched>0) 
            {
                 //then match relative positions
                 mean[0] = flatparts[0].dyx[partids.a[i]].a[4+5*j];
                 mean[1] = flatparts[0].dyx[partids.a[i]].a[5+5*j];
                 var[0] = flatparts[0].dyx[partids.a[i]].a[6+5*j];
                 var[1] = flatparts[0].dyx[partids.a[i]].a[7+5*j];
                 //find a best one in the current neighs
                 for(k=0; k<idx.length; k++)
                 {
                     if(omatched[k]==1)
                     {
                         double prob = (0.5*(dyx[2*k]- mean[0])*(dyx[2*k]- mean[0]))/(var[0]+4.5*pow(epsilon, numlayer-1)) + 0.5*((dyx[2*k+1]- mean[1])*(dyx[2*k+1]- mean[1]))/(var[1]+4.5*pow(epsilon,numlayer-1)) ;
                         //prob = exp(-prob);
                         if(maxprob<exp(-prob))
                         {
                             maxprob = exp(-prob);
                             bestk[j] = k;
                         }
                     }
                 }

                 if(maxprob > 0.5*pow(discount,numlayer-1))
                 { //it means a match

                    //get the score and its subpart locations int cy, cx, co2, radius;
                    pscore = pscore*maxprob;
                    cy = r + dyx[2*bestk[j]];
                    cx = c + dyx[2*bestk[j]+1];
                    co2 = o2list[bestk[j]];   
                    childid = model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                    //get its children from the previous layer
                    if (numlayer>2)
                    {    
                        // find a neighborhood
                        radius = model[numlayer-2].radius; 
                        pre_radius = 2;
                        if(numlayer>3)
                            pre_radius = model[numlayer-3].radius;// - radjust[numlayer-2];
                        reidx = searchneighs(idx.a[bestk[j]], pdist, ninst, radius, pre_radius, visited);
                        // printf("current id: %d, the number of neighbors: %d", idx.a[bestk[j]], reidx.length);
                        if(reidx.length>0 && !visited[idx.a[bestk[j]]])
                            // matching it
                            cfound = topdownmatch_full(pinstances, ninst, pDim, pdist, idx.a[bestk[j]], reidx, visited, model, numlayer-1, numOrient, 
                        numsub, radjust, discount,rotation, &cnode);
                        else
                            cfound = false;
                        if(!cfound)
                        {
                            if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                            if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                            if(flag!=NULL)
                            {free(flag);flag=NULL;}
                            
                             break;
                        }

                    }
                    else{
                        
                        cfound = true;
                        
                        //visited[idx.a[bestk[j]]] = true;
                        cnode.a[0] = cy;
                        cnode.a[1] = cx;
                        cnode.a[2]= model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];
                        cnode.a[3] = numlayer-1;//its child from previous layer
                        cnode.length = 1;
                        cnode.score = maxprob; 
                        // other considerations
                        // cnode.length = 0;
                        
                    }
                    
                    if(cfound){
                        flag[j] = 1;                        
                        // add children here
                        /*pnode.a[pnode.length*pDim+0] = cy;
                        pnode.a[pnode.length*pDim+1] = cx;        
                        pnode.a[pnode.length*pDim+2] = model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j];//pid
                        pnode.a[pnode.length*pDim+3] = numlayer-1;
                        pnode.length = pnode.length+1;
                        */
                        pnode = cat(pnode, cnode, pDim);
                        pnode.score = pnode.score*cnode.score;
                    }
                 
                 }
                 //add here to speed up
                 else{
                        if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                        break;
                 }

            }
            if(omatched!=NULL)
            {free(omatched);omatched = NULL;}

        }
        
        if(flag!=NULL)
        {
            //check the matching from its children
            if(afound && sum(flag, numchild)==numchild){
                found =true;
                indicator[tindex] =1;
                tscores[tindex] = pnode.score;
                partnodes[tindex] = pnode;
                /*// mask it
                visited[iinst] = false;
                // mask its children
                for(j =0; j<numchild; j++)
                    visited[idx.a[bestk[j]]] = false;
                 */
            }
            //else found = 0;
            free(flag);flag=NULL;
        }
        //else
        //    found = 0;
        if(bestk!=NULL)
        {free(bestk);bestk = NULL;} 
        
    }
    
    //free neighbors
    if(y!=NULL)
    {   free(y); y=NULL;}
    if(x!=NULL)
    {   free(x); x= NULL;}
    if(iy!=NULL)
    {   free(iy); iy= NULL;}
    if(ix!= NULL)
    {   free(ix); ix = NULL;}
    if(o2list!= NULL)
    {   free(o2list); o2list = NULL;}
    if(dyx!= NULL)
    {   free(dyx);dyx = NULL;}
    
    
    /*check how many two-subcomposition matched*/
    nummatched = 0;
    maxprob = 0;//remember the max
    int id=0; //remember id
    /* consider suppression here */
    for(i=0;i<partids.length; i++)
    {
        if(childlens[i]==2 && tscores[i]>0) //consider two children first
        {    
            nummatched = nummatched + 1;
            if(maxprob< tscores[i]) 
            {
                maxprob = tscores[i];
                id = i;
            }
        }
    }

    // only consider the 1-sub maximum score here
    if(nummatched==0)
        id= vmax(tscores, partids.length, &maxprob);
    /*copy result to output*/
    point->length =  partnodes[id].length;
    for(int jj=0; jj<pDim*point->length; jj++)
        point->a[jj] = partnodes[id].a[jj];
    point->compid = partids.a[id] +1;//compositional model id, +1 because it starts from 0
    point->score = maxprob;
    /*free some allocation*/
    if(indicator!=NULL)
    {free(indicator);indicator = NULL;}
    if(tscores!=NULL)
    {free(tscores);tscores = NULL;}
    if(partnodes!=NULL)
    {free(partnodes);partnodes = NULL;}
    if(tdlinks!=NULL)
    {   
        free(tdlinks);
        tdlinks = NULL;
    }
    if(flatparts!=NULL)
    {free(flatparts);flatparts = NULL;}
    if(childlens!=NULL)
    {free(childlens);childlens = NULL;}
    
    return found;
}





/*TOP DOWN MATCHING WITH ROTATION, deleting matched point in the next recursive step*/
bool exactmatch(double *pinstances, int ninst, int col, int *pdist, int iinst, struct array1d idx, bool *visited,
        struct Parts *model, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node *point)
{   
    // initialize the value for return
    bool found, cfound;
    found= false; cfound = false; 
    
    int o1, r,c;
    int *y,*x,*iy,*ix, *dyx,*o2list;
    //initialize the location and orientation
    r = pinstances[px2(iinst, 1, ninst, pDim)];
    c = pinstances[px2(iinst, 2, ninst, pDim)];
    o1 = pinstances[px2(iinst, 0, ninst, pDim)];
    // delete selected points by setting flag
    // pinstances[px2(iinst, 3, ninst, pDim)] =0;
    struct array1d *tdlinks = (struct array1d*)malloc(numOrient*sizeof(struct array1d));
    struct pTypes *flatparts = mappingParts(model, numOrient, numlayer, numsub, rotation, tdlinks);
    
    //get the matching parts id
    struct array1d partids = tdlinks[o1-1];
    if(partids.length<1){
        found=false;
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        return found;
    } 
    int pDim =4;
    // define variables
    int i,j,k, pid;
    double epsilon = 1.5;
    double mean[2], var[2];
    int childid, numchild, anchorid;
    int cy, cx, co2, radius, pre_radius;//for its children
    struct array1d reidx;

    /*
    // show links
    printf("---------------- begin -----------\n");
    for(i=0;i<numOrient;i++)
    {
        printf("the %d-th orient with length %d is: \n", i,tdlinks[i].length);
        for(j=0;j<tdlinks[i].length;j++)
            printf(" %d", tdlinks[i].a[j]);
        printf("end\n", i);
    }
    printf("---------------- end --------------\n");*/
    //get its neighbors
    if(idx.length>0)
    {
        y = (int*)malloc(sizeof(int)*idx.length); 
        //int x[idx.length];
        x = (int*)malloc(sizeof(int)*idx.length); 
        iy = (int*)malloc(sizeof(int)*idx.length); 
        ix = (int*)malloc(sizeof(int)*idx.length); 
        o2list = (int*)malloc(sizeof(int)*idx.length); 
        dyx = (int*)malloc(2*idx.length*sizeof(int));
        for(i=0;i< idx.length; i++)
        {

            y[i] = (int)pinstances[px2(idx.a[i], 1, ninst, pDim)];//pinstances[idx.a[i]][2];
            x[i] = (int)pinstances[px2(idx.a[i], 2, ninst, pDim)];//pinstances[idx.a[i]][3]; 
            o2list[i] = (int)pinstances[px2(idx.a[i], 0, ninst, pDim)];//pinstances[idx.a[i]][1];
            iy[i] = y[i] - r + radius;
            ix[i] = x[i] -c + radius;
            dyx[2*i]= iy[i]- radius;
            dyx[2*i+1] = ix[i]-radius; 

        }

    }
    else{
        
        if(tdlinks!=NULL)
        {   
            free(tdlinks);tdlinks = NULL;
        }
        if(flatparts!=NULL)
        {free(flatparts);flatparts = NULL;}
        
        return found;
    }
    // initialize the structure and score for each matched part id
    int *indicator = (int*)malloc(partids.length*sizeof(int));
    double *tscores = (double*)malloc(partids.length*sizeof(double));
    int *childlens = (int*)malloc(partids.length*sizeof(int));
    for(i=0;i< partids.length; i++)
    {
        indicator[i] = 0;
        tscores[i] = 0;
        childlens[i] = 1;
    }
    //define the part nodes by the number of pids
    struct Node *partnodes = (struct Node*)malloc(partids.length*sizeof(struct Node));
    // processing each part type
    int tindex =-1;
    int nummatched = 0;
    double maxprob=0;
    #pragma omp parallel for 
    for(i=0; i<partids.length; i++)
    { 
        tindex = tindex + 1;
        partnodes[i].length = 0; //initialize it
        pid = flatparts[0].dyx[partids.a[i]].a[0];//start from 1
        /*
        printf("pid=%d", pid);
        for(j=0;j<model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].length; j++)
            printf(" %f ", model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[j]);
        printf("\n");
        */
        numchild = flatparts[0].dyx[partids.a[i]].a[2];
        childlens[i] = numchild;

        //define the current node from this input point
        struct Node pnode;
        pnode.a[0] = r;
        pnode.a[1] = c;        
        pnode.a[2] = pid;//partids.a[i];
        pnode.a[3] = numlayer;
        pnode.length = 1;
        pnode.score = 1;
        
        double pscore = 1;   
        
        /*--------------------------------------------------------
         * matching for its children
         *--------------------------------------------------------
         */
        
        //indicator matching for child
        int* flag = (int*)malloc(numchild*sizeof(int));
        initialize(flag, numchild);
        int *bestk = (int*)malloc(numchild*sizeof(int));   
        
        //mask current point
        visited[iinst] = true;
        // check subparts with its neighborhood
        for(j= 0; j< numchild; j++)
        {
            
            cfound = false;
            //model[numlayer-1].ptypes_[numsub-1].dyx[pid].a[0];
            childid = flatparts[0].dyx[partids.a[i]].a[3+5*j];
            // define the children node
            struct Node cnode;
            cnode.length = 0;
            cnode.score = 1;
            // check its types
            int* omatched = find(o2list, idx.length, childid, 'e');;
            nummatched = sum(omatched, idx.length);
            //check its relative positions
            maxprob = 0;// for best match from a bunch of instances
            if(nummatched>0) 
            {
                 //then match relative positions
                 mean[0] = flatparts[0].dyx[partids.a[i]].a[4+5*j];
                 mean[1] = flatparts[0].dyx[partids.a[i]].a[5+5*j];
                 var[0] = flatparts[0].dyx[partids.a[i]].a[6+5*j];
                 var[1] = flatparts[0].dyx[partids.a[i]].a[7+5*j];
                 //find a best one in the current neighs
                 for(k=0; k<idx.length; k++)
                 {
                     if(omatched[k]==1)
                     {
                         double prob = (0.5*(dyx[2*k]- mean[0])*(dyx[2*k]- mean[0]))/(var[0]+4.5*pow(epsilon, numlayer-1)) + 0.5*((dyx[2*k+1]- mean[1])*(dyx[2*k+1]- mean[1]))/(var[1]+4.5*pow(epsilon,numlayer-1)) ;
                         //prob = exp(-prob);
                         if(maxprob<exp(-prob))
                         {
                             maxprob = exp(-prob);
                             bestk[j] = k;
                         }
                     }
                 }

                 if(maxprob > 0.5*pow(discount,numlayer-1))
                 { //it means a match

                     //get the score and its subpart locations int cy, cx, co2, radius;
                     pscore = pscore*maxprob;
                     cy = r + dyx[2*bestk[j]];
                     cx = c + dyx[2*bestk[j]+1];
                     co2 = o2list[bestk[j]];   
                     //orignal part id
                     childid = model[numlayer-1].ptypes_[numsub-1].dyx[partids.a[i]].a[3+5*j]; 
                     reidx.length =0;
                     //get its children from the previous layer
                     if (numlayer>2)
                     {    
                         // find a neighborhood
                         radius = model[numlayer-2].radius; 
                         pre_radius = 2;
                         if(numlayer>2)
                             pre_radius = model[numlayer-3].radius;// - radjust[numlayer-2];
                         reidx = searchneighs(idx.a[bestk[j]], pdist, ninst, radius, pre_radius, visited);
                         //printf("current id: %d, the number of neighbors: %d, pid=%d\n", idx.a[bestk[j]], reidx.length, childid);
                         if(reidx.length>0 && !visited[idx.a[bestk[j]]])
                             // matching it
                             cfound = hiertpdmatch(pinstances, ninst, pDim, pdist, idx.a[bestk[j]], reidx, visited, model,childid, numlayer-1, numOrient, 
                         numsub, radjust, discount,rotation, &cnode);
                     }
                    
                     if(cfound){
                        
                         flag[j] = 1;                    
                         // add children here
                         //printf("find a child\n");
                         pnode = cat(pnode, cnode, pDim);
                         pnode.score = pnode.score*cnode.score;
                     }
                     else
                     {   //add here to speed up
                          if(omatched!=NULL)
                          {free(omatched);omatched = NULL;}
                          if(bestk!=NULL)
                          {free(bestk);bestk = NULL;}
                          if(flag!=NULL)
                          {free(flag);flag=NULL;}
                            
                          break;
                     }
                     
                 }
                 else{
                     
                     if(omatched!=NULL)
                     {free(omatched);omatched = NULL;}
                     break; //for speed up
                 }
            }
            else
            {   //add here to speed up
                  if(omatched!=NULL)
                  {free(omatched);omatched = NULL;}
                  if(bestk!=NULL)
                  {free(bestk);bestk = NULL;}
                  if(flag!=NULL)
                  {free(flag);flag=NULL;}

                  break;
            }
        }
        //end of childs comparison 
        if(flag!=NULL)
        {
            //check the matching from its children
            if(sum(flag, numchild)==numchild){
                found =true;
                indicator[tindex] =1;
                tscores[tindex] = pnode.score;
                partnodes[tindex] = pnode;
                // mask it
                visited[iinst] = false;
                // mask its children
                for(j =0; j<numchild; j++)
                    visited[idx.a[bestk[j]]] = false;
            }
            //else found = 0;
            free(flag);flag=NULL;
        }
        else
            found = false;
        if(bestk!=NULL)
        {free(bestk);bestk = NULL;} 
        
    }
    
    //free neighbors
    if(y!=NULL)
    {   free(y); y=NULL;}
    if(x!=NULL)
    {   free(x); x= NULL;}
    if(iy!=NULL)
    {   free(iy); iy= NULL;}
    if(ix!= NULL)
    {   free(ix); ix = NULL;}
    if(o2list!= NULL)
    {   free(o2list); o2list = NULL;}
    if(dyx!= NULL)
    {   free(dyx);dyx = NULL;}
    
    
    /*check how many two-subcomposition matched*/
    nummatched = 0;
    maxprob = 0;//remember the max
    int id=0; //remember id
    /* consider suppression here */
    for(i=0;i<partids.length; i++)
    {
        if(childlens[i]==2 && tscores[i]>0) //consider two children first
        {    
            nummatched = nummatched + 1;
            if(maxprob< tscores[i]) 
            {
                maxprob = tscores[i];
                id = i;
            }
        }
    }

    // only consider the 1-sub maximum score here
    if(nummatched==0)
        id= vmax(tscores, partids.length, &maxprob);
    /*copy result to output*/
    point->length =  partnodes[id].length;
    for(int jj=0; jj<pDim*point->length; jj++)
        point->a[jj] = partnodes[id].a[jj];
    point->compid = partids.a[id] +1;//compositional model id, +1 because it starts from 0
    point->score = maxprob;
    /*free some allocation*/
    if(indicator!=NULL)
    {free(indicator);indicator = NULL;}
    if(tscores!=NULL)
    {free(tscores);tscores = NULL;}
    if(partnodes!=NULL)
    {free(partnodes);partnodes = NULL;}
    if(tdlinks!=NULL)
    {   
        free(tdlinks);
        tdlinks = NULL;
    }
    if(flatparts!=NULL)
    {free(flatparts);flatparts = NULL;}
    if(childlens!=NULL)
    {free(childlens);childlens = NULL;}
    
    return found;
}

/*match specific subchildren from top down hierarchy*/
bool hiertpdmatch(double *pinstances, int ninst, int col, int *pdist, int iinst, struct array1d idx, bool *visited,
        Parts *model, int pid, int numlayer, int numOrient, int numsub, double *radjust, double discount,int rotation, Node * point)
{
    
    bool found,afound, cfound;
    found= false;  
    afound = false;
    cfound = false;
    int o1, r,c, i, j, k;
    int *y,*x,*iy,*ix, *dyx,*o2list;
    // define its neighbors
    int cy, cx, co2, radius, pre_radius;//for its children
    struct array1d reidx;    
    // initialize the location and orientation
    r = pinstances[px2(iinst, 1, ninst, pDim)];
    c = pinstances[px2(iinst, 2, ninst, pDim)];
    o1 = pinstances[px2(iinst, 0, ninst, pDim)];
    radius = model[numlayer-1].radius; 
    // get its neighbors
    if(idx.length>0)
    {
        y = (int*)malloc(sizeof(int)*idx.length); 
        //int x[idx.length];
        x = (int*)malloc(sizeof(int)*idx.length); 
        iy = (int*)malloc(sizeof(int)*idx.length); 
        ix = (int*)malloc(sizeof(int)*idx.length); 
        o2list = (int*)malloc(sizeof(int)*idx.length); 
        dyx = (int*)malloc(2*idx.length*sizeof(int));
        for(i=0;i< idx.length; i++)
        {

            y[i] = (int)pinstances[px2(idx.a[i], 1, ninst, pDim)];//pinstances[idx.a[i]][2];
            x[i] = (int)pinstances[px2(idx.a[i], 2, ninst, pDim)];//pinstances[idx.a[i]][3]; 
            o2list[i] = (int)pinstances[px2(idx.a[i], 0, ninst, pDim)];//pinstances[idx.a[i]][1];
            iy[i] = y[i] - r + radius;
            ix[i] = x[i] -c + radius;
            dyx[2*i]= iy[i]- radius;
            dyx[2*i+1] = ix[i]-radius; 

        }

    }

    int pDim = 4;
    double epsilon = 1.5;
    double score = 1;
    double maxprob = 0;
    double mean[2], var[2];
    
    // printf("from model: numchild= %lf\n", model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[2]);
    int numchild = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[2];
    //printf("numlayer = %d, pid=%d, numchild= %d\n", numlayer, pid, numchild);
    
    int *flag = (int*)malloc(sizeof(int)*numchild);
    int *bestk = (int*)malloc(sizeof(int)*numchild); // for the best id in dyx vector
    int nummatched=0;
    
    //define the current node from this input point
    struct Node pnode;
    pnode.a[0] = r;
    pnode.a[1] = c;        
    pnode.a[2] = pid;//partids.a[i];
    pnode.a[3] = numlayer;
    pnode.length = 1;
    pnode.score = 1;
    pnode.compid = pid;
    double pscore = 1;   
    /*
    // add children here
    pnode->a[pnode->length*pDim+0] = r;
    pnode->a[pnode->length*pDim+1] = c;        
    pnode->a[pnode->length*pDim+2] = pid;
    // printf("r=%d, c =%d, numlayer = %d, pid=%d\n", r, c, numlayer, pid);
    pnode->a[pnode->length*pDim+3] = numlayer;
    pnode->length = pnode->length+1;
    double pscore = pnode->score;
    */
    int childid, ocid, anchodid;
    // matching anchorid
    struct Node anode;
    anode.length = 0; //initialize it
    anode.score = 1;
    anchodid = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[1];
    if (numlayer>2)
    {    
        radius = model[numlayer-2].radius; 
        pre_radius = 2;
        if(numlayer>2)
            pre_radius = model[numlayer-3].radius;// - radjust[numlayer-2];
        reidx = searchneighs(iinst, pdist, ninst, radius, pre_radius, visited);
        // printf("current id: %d, the number of neighbors: %d", idx.a[bestk[j]], reidx.length);
        if(reidx.length>0 && !visited[iinst])
            afound = hiertpdmatch(pinstances, ninst, pDim, pdist, iinst, reidx, visited, model, anchodid, numlayer-1, numOrient, 
                        numsub, radjust, discount,rotation, &anode);
    }
    else
    {
                        
        afound = true;
        /*
        //visited[idx.a[bestk[j]]] = true;
        anode.a[0] = cy;
        anode.a[1] = cx;
        anode.a[2]= childid;//model[numlayer-1].ptypes_[numsub-1].dyx[pid].a[3+5*j];
        anode.a[3] = numlayer-1;//its child from previous layer
        anode.length = 1;
        anode.score = pscore; */           
    }
    if(afound)
    {    
         cat(&pnode, &anode, pDim);
         //pnode->score = pnode->score*cnode.score;
         pnode.score = pnode.score*anode.score;
    }
    /* match two childrens */
    for(j=0; j< numchild; j = j+1)
    {   
         cfound = false;
         flag[j] = 0; //initialize the matching flag
         struct Node cnode;
         cnode.length = 0;
         cnode.score = 1;
         childid = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[3+5*j];
         ocid = id2orient(childid, model, numlayer, numsub);
         //add rotation
         ocid = ocid + rotation;
         if(ocid<1) ocid = ocid + numOrient;
         if(ocid>numOrient) ocid = ocid - numOrient;
         /*find the matched child orientation */
         int* omatched = find(o2list, idx.length, ocid, 'e');
         /*THE NUM OF MATCHED*/
         nummatched = sum(omatched, idx.length);
         maxprob = 0;// for best match for a given pid
         if(nummatched>0) 
         {
             //double* selected_dyx = selection(dyx, numdyx, indicator);
             //double* selected_o2 = selection(o2list, numdyx, indicator);
             //then match relative positions
             mean[0] = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[4+5*j];
             mean[1] = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[5+5*j];
             var[0] = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[6+5*j];
             var[1] = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[7+5*j];
             
             for(k=0; k<idx.length; k++)
             {
                 if(omatched[k]==1)
                 {
                     double prob = (0.5*(dyx[2*k]- mean[0])*(dyx[2*k]- mean[0]))/(var[0]+2.0+pow(epsilon, numlayer)) + 0.5*((dyx[2*k+1]- mean[1])*(dyx[2*k+1]- mean[1]))/(var[1]+2.0+pow(epsilon,numlayer)) ;
                     //prob = exp(-prob);
                     if(maxprob<exp(-prob))
                     {
                         maxprob = exp(-prob);
                         bestk[j] = k;
                     }
                 }
             }

             if(maxprob > 0.01*0.5*pow(discount,numlayer-1))
             { //it means a match
                 
                    //get the score and its subpart locations int cy, cx, co2, radius;
                    pscore = pscore*maxprob;
                    cy = r + dyx[2*bestk[j]];
                    cx = c + dyx[2*bestk[j]+1];
                    co2 = o2list[bestk[j]];   
                    childid = model[numlayer-1].ptypes_[numsub-1].dyx[pid-1].a[3+5*j];
                    //get its children from the previous layer
                    if (numlayer>2)
                    {    
                        // find a neighborhood
                        radius = model[numlayer-2].radius; 
                        pre_radius = 2;
                        if(numlayer>2)
                            pre_radius = model[numlayer-3].radius;// - radjust[numlayer-2];
                        reidx = searchneighs(idx.a[bestk[j]], pdist, ninst, radius, pre_radius, visited);
                        // printf("current id: %d, the number of neighbors: %d", idx.a[bestk[j]], reidx.length);
                        if(reidx.length>0 && !visited[idx.a[bestk[j]]])
                            // matching it
                            cfound = hiertpdmatch(pinstances, ninst, pDim, pdist, idx.a[bestk[j]], reidx, visited, model,childid, numlayer-1, numOrient, 
                        numsub, radjust, discount,rotation, &cnode);
                        else
                            cfound = false;
                        if(!cfound)
                        {
                            if(omatched!=NULL)
                            {free(omatched);omatched = NULL;}
                            if(bestk!=NULL)
                            {free(bestk);bestk = NULL;}
                            if(flag!=NULL)
                            {free(flag);flag=NULL;}
                            
                             break;
                        }

                    }
                    else{
                        
                        cfound = true;
                        
                        //visited[idx.a[bestk[j]]] = true;
                        cnode.a[0] = cy;
                        cnode.a[1] = cx;
                        cnode.a[2]= childid;//model[numlayer-1].ptypes_[numsub-1].dyx[pid].a[3+5*j];
                        cnode.a[3] = numlayer-1;//its child from previous layer
                        cnode.length = 1;
                        cnode.score = pscore; 
                        // other considerations
                        // cnode.length = 0;
                        
                    }
                    
                    
                    if(cfound){
                        flag[j] = 1;
                        //printf("match child, layer=%d\n", numlayer);
                        // add children here
                        cat(&pnode, &cnode, pDim);
                        //pnode->score = pnode->score*cnode.score;
                        pnode.score = pnode.score*cnode.score;
                    }
                 
             }

         }
         //no match
         else{
             
            if(omatched!=NULL)
            {free(omatched);omatched = NULL;}
            if(bestk!=NULL)
            {free(bestk);bestk = NULL;}
            if(flag!=NULL)
            {free(flag);flag=NULL;}
            break;
         }

    }//finshing of consideration of all childrens
        
    if(flag!=NULL)
    {       /*for(j=1; j<numchild; j++)
                printf("%d ",flag[j]);
            printf("\n");*/
            //check the matching from its children
            if(afound && sum(flag, numchild)>=numchild){
                found =true;
                //printf("numchild= %d, all childrens matched, numlayer=%d\n", numchild, numlayer);
                // mask it
                visited[iinst] = true;
                // mask its children
                for(j =0; j<numchild; j++)
                    visited[idx.a[bestk[j]]] = false;
            }
            else found = 0;
            free(flag);flag=NULL;
     }

     if(bestk!=NULL)
     {free(bestk);bestk = NULL;} 
    
    
    //free neighbors
    if(y!=NULL)
    {   free(y); y=NULL;}
    if(x!=NULL)
    {   free(x); x= NULL;}
    if(iy!=NULL)
    {   free(iy); iy= NULL;}
    if(ix!= NULL)
    {   free(ix); ix = NULL;}
    if(o2list!= NULL)
    {   free(o2list); o2list = NULL;}
    if(dyx!= NULL)
    {   free(dyx);dyx = NULL;}

    point->length =  pnode.length;
    for(int jj=0; jj<pDim*point->length; jj++)
        point->a[jj] = pnode.a[jj];
    point->compid = pid;//compositional model id, +1 because it starts from 0
    point->score = pnode.score;
    return found;
    
}





