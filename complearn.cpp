#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include <omp.h>  //for multi-threding
#include "complearn.h"
#include "checkparts.h"
/*define the total number of layers*/
// int totLayers = 7;
// int numsub = 2;
//#define px2(x, y, lengthx, lengthy) ((x) + (y)*(lengthx) - 0) 
//#define px3(x, y, z, lengthx, lengthy, lengthz) (z+ ((x) + ((y))*(lengthx))*lengthz - 0) 
// struct Parts model[totLayers];

void learningcomposition(struct Data *examples, int ndata, struct Parts *model,int numlayer, int numsub,  int numOrient, int radius, 
        double percent, double *param, double discount, bool flag)
{
    // learning composition
    double radjust[] = {0 , 0 , 2, 4, 8 , 10 ,15, 20};
    int neighbors[]= {-1,0, 0,-1, 1,0, 0,1};
    int ntypes;
    if(numlayer<2)
        ntypes = numOrient;
    else{
        if(flag)
            ntypes = model[numlayer-1].ptypes_[numsub].ntypes;
        else
            ntypes = model[numlayer-1].ptypes_[numsub-1].ntypes;
    }
    int numcomp = ntypes*ntypes;
    int diameter = 2*radius-1;
    int i,j,k, imidx;
    int pre_radius;
    bool compress = true;
    /// template code, may not work
    /*
    int ***sumhists;
    sumhists = (int***)malloc(diameter * sizeof(int **));
    for (i = 0; i < diameter; ++i) {
        sumhists[i] = (int**)malloc(diameter * sizeof(int *));
        for (j= 0; j < diameter; ++j)
            sumhists[i][j] = (int*)malloc(numcomp * sizeof(int));
    } 
    */
    int *sumhists = (int*)calloc(diameter * diameter * numcomp, sizeof(int));
    /*
    //initialize it
    #pragma omp parallel for private(i,j,k), shared(sumhists)
    for(i=0;i<diameter; i++)
        for(j=0;j<diameter; j++)
        for(k=0;k<numcomp; k++){
            sumhists[px3(i,j,k,diameter, diameter, numcomp)] = 0;
            //sumhists[i][j][k] = 0;
        }
    */
    // copy data

    // for the time computing
    clock_t tstart, tend; 
    double elapsed_time;
    tstart = clock() ;
    
    
    int r, c, o1, o2;
    struct array1d idx, reidx;
    /*=========================================================
     * accumulate histogram
     *=========================================================
     */
    int **pdist = (int**)malloc(ndata*sizeof(int*));
    // double**pdist = (double**)malloc(ndata*sizeof(double*));
    int* prob = (int*)calloc(numcomp, sizeof(int)); // allocate space and initialize it as zeros
    int* morients;
    int imh,imw,col,ninst;
    double *pinstances;
    //#pragma omp parallel for private(imidx) schedule(dynamic) reduction(+:sumhists)
    //#pragma omp parallel //for private(imidx, morients, imh,imw,col,ninst), shared(pdist, sumhists)
    for(imidx=0; imidx< ndata; imidx++)
    {
        pdist[imidx] = NULL;
        morients = (int*)examples[imidx].morients;
        imh = examples[imidx].imh;
        imw = examples[imidx].imw;
        ninst = examples[imidx].pinstances[numlayer-1].row;
        if(ninst<1)
        { 
            continue;
        }
        col = examples[imidx].pinstances[numlayer-1].col;
        pinstances = examples[imidx].pinstances[numlayer-1].elements;
        //printf("image number: %d, instance number: %d, col: %d\n", imidx, ninst, col);
        /*compute the distance between points*/
        if(numlayer>3)
        {
            pdist[imidx] = (int*)malloc(ninst*ninst*sizeof(int));
            funDist(pinstances, ninst, col, pdist[imidx]);
        }
        else{
            
            pdist[imidx] = (int*)calloc(imh*imw*MDIM, sizeof(int));
            //initmask(pdist[imidx], imh, imw);
            buildmask3(pdist[imidx],imh, imw, MDIM, pinstances, ninst, col);
        }
        /*compute the neighborhood for each point*/
        int iinst =0; 
        int y, x, iy, ix, olist;
        //#pragma omp parallel //for private(iinst,y,x,iy,ix,olist), shared(pdist)
        for(iinst = 0; iinst< ninst; iinst++)
        {

            r = pinstances[px2(iinst, 1, ninst, col)];
            c = pinstances[px2(iinst, 2, ninst, col)];
            o1 = pinstances[px2(iinst, 0, ninst, col)];
            //if(numlayer>1)
            //struct cellData ele = examples[imidx].lambda[numlayer-1].elements[iinst];
            // find the neighbors 
            //idx1 = searchneighs(r, c, pinstances, ninst, col, radius, numlayer);
            pre_radius = 4;
            if(numlayer>1)
                pre_radius = model[numlayer-1].radius - radjust[numlayer-1];
            //idx2 = searchneighs2(r, c, pinstances, ninst, col, pre_radius, numlayer-1, radjust);
            //idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer);  
            //idx = setdiff(idx1,idx2);
            if(numlayer>3)
                idx = searchneighs(iinst, pdist[imidx], ninst, radius, pre_radius, numlayer);
            else
                idx = lookneighs3(r,c, pdist[imidx], imh, imw, MDIM, radius, pre_radius);
            if(idx.length> 0)
            {   
                for(i=0;i< idx.length; i++)
                {
                    y = pinstances[px2(idx.a[i], 1, ninst, col)];//pinstances[idx.a[i]][2];
                    x = pinstances[px2(idx.a[i], 2, ninst, col)];//pinstances[idx.a[i]][3]; 
                    o2 = pinstances[px2(idx.a[i], 0, ninst, col)];//pinstances[idx.a[i]][1];
                    iy = y - r + radius;
                    ix = x - c + radius;
                    
                    #pragma omp critical
                    sumhists[px3(iy-1, ix-1, (o1-1)*ntypes+o2-1, diameter,diameter, numcomp)]= sumhists[px3(iy-1, ix-1, (o1-1)*ntypes+o2-1, diameter, diameter, numcomp)] + 1;
                    prob[(o1-1)*ntypes+o2-1] += 1;     
                }
                
            }

        }
        //free the allocation for distance computation
        //free(pdist);

    }

    
    elapsed_time = (tend-tstart)/(double)CLOCKS_PER_SEC;
    printf("the time computing for histogram: %f\n", elapsed_time);
    // count # of learned parts 
    int nparts = 0;
    int pDim = 9;
    double *ptypes = (double*)malloc(pDim*MAX_NUM_PARTS*sizeof(double));
    int* nptypes = (int*)malloc(numsub*sizeof(int));
    double* hist = (double*)malloc(diameter*diameter*sizeof(double));
    int anchorid, cid; 
    bool mflag = true;
    
    char strname[45];
    
    // find the significant compositions 
    int np = 2;  //how many parts for each histogram
    int numPoints = 0; //counting how many points in each center
    double *dataPoints;
    double *final;
    int numpeak = 0; //counting the total number of peaks
    //#pragma omp parallel for private(i,dataPoints,numPoints,final,numpeak) schedule(dynamic) reduction(+:nparts), shared(ptypes)
    for(i=0; i< numcomp; i++)
    {   
         //consider speed up here
         if(prob[i]< param[numlayer-1]*ndata)
             continue;
         
         if(nparts > MAX_NUM_PARTS+EXT_NUM_PARTS-1-2*np)
             break;
         
         type2orient(i+1, ntypes, &anchorid, &cid);
         // initialize the histogram
         #pragma omp parallel for private(j,k)
         for(j= 0; j<diameter; j++)
            for(k=0; k<diameter; k++)
            {
                hist[j*diameter+k] = 0;
            }
         
         #pragma omp parallel for private(j,k),shared(i,sumhists)
         for(j= 0; j<diameter; j++)
            for(k=0; k<diameter; k++)
            {    //hist[j*diameter+k] 
                  hist[px2(j, k, diameter, diameter)] =  sumhists[px3(j, k, i, diameter, diameter, numcomp)];
                //hist[j*diameter+k]  =  sumhists[j][k][i]; 
         }
         
         //sprintf(strname, "%s_%d_%d_%d_", "hist", numlayer, anchorid, cid,".txt");
         //writeArraytoFile(hist, diameter, diameter,strname);
         //prob[i] = sum(hist, diameter*diameter);
         
         dataPoints = NULL;     
         // int numPoints = 0;
         dataPoints = hist2points(hist, diameter, &numPoints);
         /*draw the histogram here*/
         /*if(numPoints>0){
            sprintf(strname, "%s_%d_%d_%d_%s", "hist", numlayer, anchorid, cid,".txt");
            writeArraytoFile(hist, diameter, diameter,strname);
         }*/
         // int numpeak = 0;
         final = local_maxima(hist,radius, &numpeak, neighbors);
         /*
         if(anchorid==1&&cid ==24)
         for(j=0;j<numpeak; j++){
             printf("the %d-th value: %f, coordinate: %f, %f\n",j,final[j*3+0],final[j*3+1], final[j*3+2]);
             if(numPoints>0){
                sprintf(strname, "%s_%d_%d_%d_%s", "hist", numlayer, anchorid, cid,".txt");
                writeArraytoFile(hist, diameter, diameter,strname);
             }
         }*/
         //int np = sizeof(final)/(sizeof(double)*3);
         //copy the selected points
         int* subdataPoints=NULL;
         double mean[2],var[2];
         //#pragma omp parallel for private(j) schedule(dynamic) reduction(+:nparts), shared(final, dataPoints)
         for(j= 0; j< MIN(numpeak,np); j++)
         { 
            r = final[j*3+1];
            c = final[j*3+2];
            struct cellData cidx = searchneighs_ex(r, c, dataPoints, prob[i], 2, 0.4*radius, numlayer);
            //subdataPoints = getneighs(r, c, hist_bk, diameter, diameter, round(0.4*radius), &numPoints);
            // printf("the current points #: %d\n", numPoints);
            if(cidx.length> param[numlayer-1]*ndata)
            //if(final[j*3+0]>param[numlayer])
            {
                
                mean[0] = 0; mean[1] = 0;
                var[0] = 0; var[1] = 0;
                //double* mean = (double*)malloc(2*sizeof(double));
                //double* var = (double*)malloc(2*sizeof(double));
                subdataPoints = NULL;
                subdataPoints = subcopy(dataPoints, numPoints, cidx, 2);
                //compute mean and variance
                //meanvariance(subdataPoints, numPoints, mean, var);
                meanvariance(subdataPoints, cidx.length, mean, var, numlayer);
                //int dyx[8];
                if(nparts>MAX_NUM_PARTS-1 && mflag)
                {
                    // only realloc the memory once, mflag controls it
                    double* extptypes = (double*)realloc(ptypes, pDim*(MAX_NUM_PARTS+EXT_NUM_PARTS)*sizeof(double));
                    if (extptypes!=NULL) 
                    {
                        ptypes=extptypes;
                        mflag = false;

                    }
                    else
                    {
                        printf("Allocation error\n");
                        exit(0);
                    }
            
                }

                ptypes[nparts*pDim+0] = i;
                ptypes[nparts*pDim+1] = anchorid;
                ptypes[nparts*pDim+2] = 1; 
                ptypes[nparts*pDim+3] = cid;
                ptypes[nparts*pDim+4] = mean[0]-radius;
                ptypes[nparts*pDim+5] = mean[1]-radius;
                ptypes[nparts*pDim+6] = var[0];
                ptypes[nparts*pDim+7] = var[1];
                ptypes[nparts*pDim+8] = cidx.length;
                // increase the number of parts
                nparts++;   
                if(subdataPoints)
                {free(subdataPoints);subdataPoints=NULL;}
            }
            //free it after use it
            free(cidx.a); 
            
         }
         /*if(nparts > MAX_NUM_PARTS+EXT_NUM_PARTS-1-2*np)
             break;*/
         if(dataPoints) {free(dataPoints);dataPoints = NULL;}
         if(final) {free(final);final = NULL;}

    }/*statistical parts*/
    printf("the number of 1-sub parts learned: %d\n", nparts);
    
    
    if(hist)
    {
        free(hist);
        hist=NULL;
    }
    if(prob)
    {
        free(prob);
        prob = NULL;
    }
    
    /*
    for (i = 0; i < diameter; ++i)
        for (j= 0; j < diameter; ++j)
            free(sumhists[i][j]);
    for (i = 0; i < diameter; ++i)    
        free(sumhists[i]);
    */
    if(sumhists){
        free(sumhists);
        sumhists = NULL;
    }
/*=======================================================
 * store one composition in the model
 *=======================================================
 */
    
    nptypes[0] = nparts;//store all the 1-sub parts learned
    
    
    double sumpprob = 0;
    double threshold =0;
    
    double* pprob = (double*)malloc(nparts*sizeof(double));
    /*set suppression weight = 0.4*/
    double *weights = setweights(ptypes, nptypes, pDim, 1, 0.4);
    #pragma omp parallel for private(i) schedule(dynamic) reduction(+:sumpprob)
    for(i=0; i< nptypes[0]; i++)
    {   //adjust weights here
        //ptypes[i*pDim+8] = ptypes[i*pDim+8]*weights[i];
        pprob[i] = ptypes[i*pDim+8];
        sumpprob = sumpprob + pprob[i];
    }
    #pragma omp parallel for
    for(i=0; i< nptypes[0]; i++)
    {   
        pprob[i]= -pprob[i]/(0.1+sumpprob);
	
    }
    //quickSort2(pprob,0, nptypes[0]-1);
    bubblesort(pprob, nptypes[0]);
    //compute the threshold for part selection
    if(nptypes[0] > MAX_NUM_PARTS-1){
        
        threshold = -pprob[MAX_NUM_PARTS-1]*sumpprob;
    }
    else{
        threshold = 0; 
    }
    //start to select higher freqency parts
    nparts = 0; 
    sumpprob = 0;
    //#pragma omp parallel for private(i) schedule(dynamic) reduction(+:nparts)
    for(i=0; i<nptypes[0]; i++) //only copy part of the model 
    {
            
        if(ptypes[i*pDim+8]>threshold)
        {     
            model[numlayer].ptypes_[0].dyx[nparts].a[0] = nparts + 1; //i+1;//ptypes[i*pDim+0]+1;// id+1 because it starts from 0
            model[numlayer].ptypes_[0].dyx[nparts].a[1] = ptypes[i*pDim+1]; // anchor id
            model[numlayer].ptypes_[0].dyx[nparts].a[2] = ptypes[i*pDim+2]; // # of child
            model[numlayer].ptypes_[0].dyx[nparts].a[3] = ptypes[i*pDim+3]; // child id
            model[numlayer].ptypes_[0].dyx[nparts].a[4] = ptypes[i*pDim+4];
            model[numlayer].ptypes_[0].dyx[nparts].a[5] = ptypes[i*pDim+5];
            model[numlayer].ptypes_[0].dyx[nparts].a[6] = ptypes[i*pDim+6];
            model[numlayer].ptypes_[0].dyx[nparts].a[7] = ptypes[i*pDim+7];
            model[numlayer].ptypes_[0].dyx[nparts].a[8] = ptypes[i*pDim+8];      
            model[numlayer].ptypes_[0].dyx[nparts].length = 8;   
            /*remember the probability*/
            sumpprob = ptypes[i*pDim+8] + sumpprob;
            nparts = nparts + 1;
        }
            
        if(nparts>MAX_NUM_PARTS-1)
              break;
        
    }
    printf("select one subpart = %d, threshold %f\n", nparts, threshold);
    // relabel its ID and normalize pprob
    //for(i = 0; i<nptypes[0];i++)
    #pragma omp parallel for
    for(i = 0; i<nparts;i++)   
    {  
        //model[numlayer].ptypes_[0].dyx[i].a[0] = i+1; //renew its id because it might larger than MAX_NUM_PARTS
        model[numlayer].ptypes_[0].p.a[i]= double(model[numlayer].ptypes_[0].dyx[i].a[8])/(sumpprob+0.01);
    }
    model[numlayer].ptypes_[0].p.length = nparts;
    model[numlayer].ptypes_[0].ntypes = nparts;
    printf("the number of probability 1-sub parts learned: %d\n", nparts);
    struct array1d *lutlinks = NULL;
    
    lutlinks =  (struct array1d *)malloc(ntypes*sizeof(array1d));        
    //initialize look up tables
    for(i = 0; i< ntypes; i++)
    {    
        // anchorid =(int)model[numlayer-1].ptypes_[numsub-1].dyx[i].a[1]; 
        lutlinks[i].length = 0;
    }
    if(compress)
    {    
        nptypes[0] = nparts;
        //build links/look up tables to speed up
        int clen = 0;
        // compute links from sanja's paper // 
        for(i=0; i< model[numlayer].ptypes_[0].ntypes; i++) //process current part types
        {
            anchorid = (int)model[numlayer].ptypes_[0].dyx[i].a[1];
            clen = lutlinks[anchorid-1].length;
            if(clen>MAX_LEN-2)
                continue;
            lutlinks[anchorid-1].a[clen] = model[numlayer].ptypes_[0].dyx[i].a[0]-1;//make it start from 0
            clen = clen + 1;
            lutlinks[anchorid-1].length = clen;//model[numlayer-1].plinks_[anchorid-1].length +1;

        }
 
    }
    else
    {
        int clen = 0;
        // compute links from sanja's paper // 
        for(i=0; i< nptypes[0]; i++) //process current part types
        {
            anchorid = (int)ptypes[i*pDim +1];
            clen = lutlinks[anchorid-1].length;
            if(clen>MAX_LEN-2)
                continue;
            lutlinks[anchorid-1].a[clen] = i;
            clen = clen + 1;
            lutlinks[anchorid-1].length = clen;//model[numlayer-1].plinks_[anchorid-1].length +1;

        }
        
        
    }
    /*==================================================================
     * consider more than one-subcomposition
     *==================================================================
     */   
    
    //update the number of parts according to frequency selection
    if(compress) nptypes[0] = nparts;
    //two-subcomposition
    int subcomb = nptypes[0]*nptypes[0] + nptypes[0];
    double* lists = (double*)malloc(subcomb*sizeof(double));
    if(numsub>1)
    {
         /*initialize it*/
        #pragma omp parallel for
         for(i= 0; i< subcomb; i++)
             lists[i] = 0;
         /*inference first*/ 
         //#pragma omp parallel for private(imidx) schedule(dynamic) reduction(+:lists)
         //#pragma omp for
         for(imidx=0; imidx< ndata; imidx++)
         {
              // copy data
             /* 
             int morients[imh][imw];
              
              for(i = 0; i < imh; i++)
                    for(j=0; j< imw; j++)
                    {
                        morients[i][j] = data(imidx).morients[i][j];
                    }
              ninst = 
              int pinstances[ninst][4];
              
              int r, c, o1, o2;
              struct array1d idx, idx2;
              */
             
            morients = (int*)examples[imidx].morients;
            imh = examples[imidx].imh;
            imw = examples[imidx].imw;
            ninst = examples[imidx].pinstances[numlayer-1].row;
            if(ninst<1)
                continue;
            col = examples[imidx].pinstances[numlayer-1].col;
            pinstances = examples[imidx].pinstances[numlayer-1].elements;  
            int iinst;
            int *y, *x, *iy, *ix, *o2list, *dyx;
            #pragma omp parallel private(iinst,y,x,iy,ix,o2list,dyx), shared(pinstances,lists)
            for(iinst = 0; iinst< ninst; iinst++)
            {

                 r = pinstances[px2(iinst, 1, ninst, col)];
                 c = pinstances[px2(iinst, 2, ninst, col)];
                 o1 = pinstances[px2(iinst, 0, ninst, col)];
                 // find the neighbors 
                 //idx1 = searchneighs(r, c, pinstances, ninst, col, radius, numlayer);
                 pre_radius = 2;
                 if(numlayer>1)
                     pre_radius = model[numlayer-1].radius;// - radjust[numlayer-1];
                 //idx2 = searchneighs2(r, c, pinstances, ninst, col, pre_radius, numlayer-1, radjust);
                 //idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer);  
                 //idx = setdiff(idx1,idx2);
                 if(numlayer>3)
                     idx = searchneighs(iinst, pdist[imidx], ninst, radius, pre_radius, numlayer);
                 else
                     idx = lookneighs3(r,c, pdist[imidx],imh, imw, MDIM, radius, pre_radius);
                 if(idx.length>0)
                 {
                      y = (int*)malloc(sizeof(int)*idx.length); 
                      //int x[idx.length];
                      x = (int*)malloc(sizeof(int)*idx.length); 
                      //int iy[idx.length]; 
                      //int ix[idx.length];
                      //int olist[idx.length];
                      iy = (int*)malloc(sizeof(int)*idx.length); 
                      ix = (int*)malloc(sizeof(int)*idx.length); 
                      o2list = (int*)malloc(sizeof(int)*idx.length); 
                      dyx = (int*)malloc(2*idx.length*sizeof(int));
                      for(i=0;i< idx.length; i++)
                      {

                            y[i] = pinstances[px2(idx.a[i], 1, ninst, col)];//pinstances[idx.a[i]][2];
                            x[i] = pinstances[px2(idx.a[i], 2, ninst, col)];//pinstances[idx.a[i]][3]; 
                            o2list[i] = pinstances[px2(idx.a[i], 0, ninst, col)];//pinstances[idx.a[i]][1];
                            iy[i] = y[i] - r + radius;
                            ix[i] = x[i] -c + radius;
                            dyx[2*i]= iy[i]- radius;
                            dyx[2*i+1] = ix[i]-radius; 
                      }
                            
                      int numdyx = idx.length;
                      struct array1d ctypeids;
                      if(compress)
                            ctypeids = checksubparts(o1, r,c, dyx, o2list,numdyx, model[numlayer].ptypes_, nptypes, numlayer,numsub-1, discount, lutlinks);
                      else
                            ctypeids = checksubparts(o1, r,c, dyx, o2list,numdyx, ptypes, nptypes, pDim, numlayer,numsub-1, discount, lutlinks);
                      #pragma omp critical
                      if(ctypeids.length>1)
                      {

                           for(j=0; j< ctypeids.length-1; j++)
                               for(k=j+1; k<ctypeids.length; k++)
                               {
                                     
                                    lists[nptypes[0] + (ctypeids.a[j])*nptypes[0]+ctypeids.a[k]] = 
                                    lists[nptypes[0] + (ctypeids.a[j])*nptypes[0]+ctypeids.a[k]] +1;

                               }
                      }
                      else
                      {
                           if(ctypeids.length==1)
                                lists[ctypeids.a[0]] = lists[ctypeids.a[0]] +1;

                      }
                              

                      if(y) {free(y); y =NULL;}
                      if(x) {free(x); x = NULL;}
                      if(iy) {free(iy);iy = NULL;}
                      if(ix) {free(ix); ix = NULL;}
                      if(o2list) {free(o2list); o2list=NULL;}
                      if(dyx) {free(dyx); dyx = NULL;}
                 }
            }
         
         }

    }//end of subcomposition
    // percent = 0.95;
    /*find the significant combination from the list*/
    double* lists_bk = (double*)malloc(subcomb*sizeof(double));
    //no weight consider here
    
    double sumlists = sum(lists, subcomb);
    //printf("sum of the lists before sort: %2.1f\n", sumlists);
    #pragma omp parallel for
    for(i=0; i< subcomb; i++)
    {
        lists_bk[i] = lists[i]; //for back up
        //printf("%2f ", lists[i]);
        //if((i+1)%20==0) printf("\n");
        lists[i] = -lists[i]/sumlists; //negative for sorting
    }/*
    // consider weight here
    for(i=0; i< nptypes[0]; i++)
    {
        lists[i] = lists[i]*weights[i]; //set weight here
        lists_bk[i] = lists[i];
    }
    
    double temp = 1;
    for(i=nptypes[0]; i< subcomb; i++)
    {
        if(lists[i]<1)
        {
            lists_bk[i] = lists[i];
            continue;
        }
        type2orient(i+1-nptypes[0], nptypes[0], &anchorid, &cid);
        anchorid = anchorid -1; //start from 0;
        cid = cid -1;
        if(compress)
            temp = adjweight(model[numlayer].ptypes_, nptypes, 1, anchorid, cid, 0.4);
        else
            temp = adjweight(ptypes, nptypes, pDim, 1, anchorid, cid, 0.4);
        
        lists[i] = lists[i]*temp; //set weight here
        lists_bk[i] = lists[i];
    }
    double sumlists = sum(lists, subcomb);
    for(i=0; i< subcomb; i++){
        
        lists[i] = -lists[i]/sumlists; //negative for sorting
    }
    */
    /*quick sorting */
    //printf("lists' address before sort: %x\n", lists);
    quickSort2(lists, 0, subcomb-1);
    // bubblesort(lists, subcomb);
    //printf("lists' address after sort: %x\n", lists);
    /*accumulation */
    double* accumlists = (double*)malloc(subcomb*sizeof(double));
    
    i=0; /*initialize it first*/
    accumlists[i] = lists[0];
    for(i=1; i< subcomb; i++)
        accumlists[i] = accumlists[i-1]+lists[i];
   
    for(i=0; i< subcomb; i++)
    {
        if(-accumlists[i] > percent)
            break;
    }
    //printf("the accumulation is: %f\n", accumlists[i]);
    /*generate a threshold according to the percent*/
    if(i>subcomb-1)
        i = subcomb-1;
    threshold = -lists[i]*sumlists-0.001;
    int numof2parts = 0;
    for(i =0; i< subcomb; i++)
    {
        if(lists_bk[i] > threshold)
            numof2parts++;
    }
    nptypes[numsub-1] = numof2parts;
    printf("the number of 2-sub parts with not maximum control is: %d\n", numof2parts);
    // find the most frequent parts if exceeding define MAX number of Parts
    /*if(numof2parts >MAX_NUM_PARTS+EXT_NUM_PARTS-1)
    {
        threshold = -lists[MAX_NUM_PARTS+EXT_NUM_PARTS-1]*sumlists;
        nptypes[numsub-1] = MAX_NUM_PARTS+EXT_NUM_PARTS;
    }*/
    if(numof2parts > MAX_NUM_PARTS-1)
    {
        threshold = -lists[MAX_NUM_PARTS-1]*sumlists;
        nptypes[numsub-1] = MAX_NUM_PARTS;
    }
    
    //printf("the 2-sub threshold is: %f\n", threshold);
    
    /*save the selected composition*/
    /*one subcomposition*/
    int ncount = 0;
    if(compress)
    {   
        #pragma omp parallel for private(i) schedule(dynamic) reduction(+:ncount)
        for(i=0; i<nptypes[0]; i++) 
        {
            if(lists_bk[i] > threshold)
            {

                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[0] = model[numlayer].ptypes_[0].dyx[i].a[0];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[1] = model[numlayer].ptypes_[0].dyx[i].a[1];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[2] = model[numlayer].ptypes_[0].dyx[i].a[2]; 
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[3] = model[numlayer].ptypes_[0].dyx[i].a[3];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[4] = model[numlayer].ptypes_[0].dyx[i].a[4];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[5] = model[numlayer].ptypes_[0].dyx[i].a[5];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[6] = model[numlayer].ptypes_[0].dyx[i].a[6];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[7] = model[numlayer].ptypes_[0].dyx[i].a[7];
                // model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = ptypes[i*pDim+8];      
                model[numlayer].ptypes_[numsub-1].dyx[ncount].length = 8;   
                /*remember the probability*/
                model[numlayer].ptypes_[numsub-1].p.a[ncount]= lists_bk[i]/(sumlists+0.01);
                ncount = ncount + 1;

            }
        }
        /*two subcomposition*/
        printf("one-two subparts split point: %d, with total 1-subparts: %d\n", ncount, i);
        #pragma omp parallel for private(i) schedule(dynamic) reduction(+:ncount)
        for(i=nptypes[0]; i< subcomb; i++) 
        {
            if(lists_bk[i] > threshold)
            {

                //int anchorid, cid;
                type2orient(i+1-nptypes[0], nptypes[0], &anchorid, &cid);
                anchorid = anchorid -1; //start from 0;
                cid = cid -1;
                //printf("two subpart'ids %d and %d, index= %d\n", anchorid, cid, i);
                //printf("compid: %d, freq= %f, anchorid= %d, childid = %d\n", i,lists_bk[i], anchorid, cid);
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[0] = model[numlayer].ptypes_[0].dyx[anchorid].a[0];//pid
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[1] = model[numlayer].ptypes_[0].dyx[anchorid].a[1];//anchorid
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[2] = 2; // number of child 
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[3] = model[numlayer].ptypes_[0].dyx[anchorid].a[3];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[4] = model[numlayer].ptypes_[0].dyx[anchorid].a[4];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[5] = model[numlayer].ptypes_[0].dyx[anchorid].a[5];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[6] = model[numlayer].ptypes_[0].dyx[anchorid].a[6];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[7] = model[numlayer].ptypes_[0].dyx[anchorid].a[7];
                // model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = ptypes[anchorid*pDim+8];   

                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = model[numlayer].ptypes_[0].dyx[cid].a[3];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[9] = model[numlayer].ptypes_[0].dyx[cid].a[4];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[10] = model[numlayer].ptypes_[0].dyx[cid].a[5];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[11] = model[numlayer].ptypes_[0].dyx[cid].a[6];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[12] = model[numlayer].ptypes_[0].dyx[cid].a[7];
                //model[numlayer].ptypes_[numsub-1].dyx[ncount].a[2] = ptypes[cid*pDim+8];   



                model[numlayer].ptypes_[numsub-1].dyx[ncount].length = 13;   
                /*remember the probability*/
                model[numlayer].ptypes_[numsub-1].p.a[ncount]= lists_bk[i]/(sumlists+0.01);
                ncount = ncount + 1;

                /*if(ncount > MAX_NUM_PARTS+EXT_NUM_PARTS-1)
                     break;*/
            }
        }
    }
    else{
        for(i=0; i<nptypes[0]; i++) 
        {
            if(lists_bk[i] > threshold)
            {

                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[0] = ptypes[i*pDim+0];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[1] = ptypes[i*pDim+1];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[2] = ptypes[i*pDim+2]; 
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[3] = ptypes[i*pDim+3];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[4] = ptypes[i*pDim+4];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[5] = ptypes[i*pDim+5];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[6] = ptypes[i*pDim+6];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[7] = ptypes[i*pDim+7];
                // model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = ptypes[i*pDim+8];      
                model[numlayer].ptypes_[numsub-1].dyx[ncount].length = 8;   
                /*remember the probability*/
                model[numlayer].ptypes_[numsub-1].p.a[ncount]= lists_bk[i]/(sumlists+0.01);
                ncount = ncount + 1;

            }
        }
        /*two subcomposition*/
        printf("one-two subparts split point: %d\n", i);
        for(i=nptypes[0]; i< subcomb; i++) 
        {
            if(lists_bk[i] > threshold)
            {

                //int anchorid, cid;
                type2orient(i+1-nptypes[0], nptypes[0], &anchorid, &cid);
                anchorid = anchorid -1; //start from 0;
                cid = cid -1;
                //printf("compid: %d, freq= %f, anchorid= %d, childid = %d\n", i,lists_bk[i], anchorid, cid);
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[0] = ptypes[anchorid*pDim+0];//pid
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[1] = ptypes[anchorid*pDim+1];//anchorid
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[2] = 2; // number of child 
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[3] = ptypes[anchorid*pDim+3];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[4] = ptypes[anchorid*pDim+4];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[5] = ptypes[anchorid*pDim+5];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[6] = ptypes[anchorid*pDim+6];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[7] = ptypes[anchorid*pDim+7];
                // model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = ptypes[anchorid*pDim+8];   

                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[8] = ptypes[cid*pDim+3];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[9] = ptypes[cid*pDim+4];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[10] = ptypes[cid*pDim+5];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[11] = ptypes[cid*pDim+6];
                model[numlayer].ptypes_[numsub-1].dyx[ncount].a[12] = ptypes[cid*pDim+7];
                //model[numlayer].ptypes_[numsub-1].dyx[ncount].a[2] = ptypes[cid*pDim+8];   



                model[numlayer].ptypes_[numsub-1].dyx[ncount].length = 13;   
                /*remember the probability*/
                model[numlayer].ptypes_[numsub-1].p.a[ncount]= lists_bk[i]/(sumlists+0.01);
                ncount = ncount + 1;

                if(ncount > MAX_NUM_PARTS+EXT_NUM_PARTS-1)
                     break;
            }
        }
    
    }
    //model[numlayer].ptypes_[0].ntypes =nparts;// nptypes[0];
    model[numlayer].ptypes_[numsub-1].ntypes = ncount;//remember how many types learned
    printf("the number of 2-sub parts learned: %d\n", ncount);
    model[numlayer].ptypes_[numsub-1].p.length = ncount;
    /*sort the part from high to low*/
    
    
    /*rellocate id for each parts*/
    for(i=0; i< model[numlayer].ptypes_[numsub-1].ntypes; i++)
    {
        model[numlayer].ptypes_[numsub-1].dyx[i].a[0] = i+1; // +1 because it starts from 0
        
    }
    
    if(numlayer>0)
    {     
        
        model[numlayer-1].plinks_ =  (struct array1d *)malloc(ntypes*sizeof(array1d));
        //initialize links for previous layer
        for(i = 0; i< ntypes; i++)
        {    
            // anchorid =(int)model[numlayer-1].ptypes_[numsub-1].dyx[i].a[1]; 
            model[numlayer-1].plinks_[i].length = 0;
        }
        int clen = 0;
        /* compute links from sanja's paper */ 
        for(i=0; i< model[numlayer].ptypes_[numsub-1].ntypes; i++) //process current part types
        {
            anchorid = (int)model[numlayer].ptypes_[numsub-1].dyx[i].a[1];
            clen = model[numlayer-1].plinks_[anchorid-1].length;
            if(clen>MAX_LEN-2)
                continue;
            model[numlayer-1].plinks_[anchorid-1].a[clen] = model[numlayer].ptypes_[numsub-1].dyx[i].a[0];
            clen = clen + 1;
            model[numlayer-1].plinks_[anchorid-1].length = clen;//model[numlayer-1].plinks_[anchorid-1].length +1;
            
        }
    
    }
    model[numlayer].radius = radius;
    /*short cut*/
    for(i=0;i<numsub;i++)
    {    
        model[numlayer].ntypes_[i]= model[numlayer].ptypes_[i].ntypes;
        //printf("the layer number:%d, its corresponding part types: %d\n", numlayer, model[numlayer].ntypes_[i]);
    }
    
    if(accumlists)
    {
        free(accumlists);
        accumlists = NULL;
    }
    if(lists_bk){
        free(lists_bk);
        lists_bk = NULL;
    }
    if(lists){
        free(lists);
        lists = NULL;
    }
    if(ptypes)
    {
        free(ptypes);
        ptypes = NULL;
    }
    if(pprob){
        free(pprob);
        pprob = NULL;
    }

    if(weights)
    {
        free(weights);
        weights = NULL;
    }
    
    //free look up tables
    if(lutlinks!=NULL)
    {
        free(lutlinks);
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
    
}//endof learning composition
 

// initialize the mask
void initmask(int *img, int imh, int imw)
{
    
    for(int i=0; i< imh; i++)
        for(int j=0; j<imw;j++)
            img[px2(i,j,imh,imw)] = 0;
    
}

// initialize the mask
void initmask3(int *img, int imh, int imw, int pDim)
{
    
    for(int i=0; i< imh; i++)
        for(int j=0; j<imw;j++)
            for(int k=0; k<pDim; k++)
                img[px3(i,j,k, imh,imw, pDim)] = 0;
    
}

// an image to store all the IDs
void buildmask(int *img, int imh, int imw, double *pinstances, int ninst, int col)
{
   int i,j;
   int y,x;
   for(i=0;i< ninst; i++)
   {
       y = pinstances[px2(i,1,ninst, col)];
       x = pinstances[px2(i,2,ninst, col)];
       img[px2(y,x, imh,imw)]= i+1;
   }
           
}


// an image to store all the IDs
void buildmask3(int *img, int imh, int imw, int pDim, double *pinstances, int ninst, int col)
{
   int i,j, k;
   int y,x;
   for(i=0;i< ninst; i++)
   {
       y = pinstances[px2(i,1,ninst, col)];
       x = pinstances[px2(i,2,ninst, col)];
       k = 0;
       while(k< pDim-1 && img[px3(y,x,k, imh,imw, pDim)]>0)
           k++;
       
       img[px3(y, x, k, imh, imw, pDim)] = i+1;
       
           
   }
           
}

struct array1d lookneighs(int y, int x, int *img, int imh, int imw, int radius, int pre_radius)
{ 
    
   int i,j;
   struct array1d idx;
   idx.length = -1;
   //crop the windows to get the neighbors id
   int miny = MAX(y-radius+1, 1);
   int minx = MAX(x-radius+1, 1);
   int maxy = MIN(y+radius-1, imh);
   int maxx = MIN(x+radius-1, imw);

   int pre_miny = MAX(y-pre_radius+1, 1);
   int pre_minx = MAX(x-pre_radius+1, 1);
   int pre_maxy = MIN(y+pre_radius-1, imh);
   int pre_maxx = MIN(x+pre_radius-1, imw);   
   
   for(i= miny; i<= maxy; i++)
       for(j= minx; j<= maxx; j++)
       {
           
            if(i>=pre_miny && i<=pre_maxy && j>=pre_minx && j<=pre_maxx)
                continue;
            if(img[px2(i,j, imh,imw)]>0){
                
               idx.length = idx.length + 1;
               idx.a[idx.length] = img[px2(i,j, imh,ihw)]-1;
    
           }
       }
   idx.length = idx.length + 1;
   return idx;
}


struct array1d lookneighs3(int y, int x, int *img, int imh, int imw, int pDim, int radius, int pre_radius)
{ 
    
   int i,j,k;
   struct array1d idx;
   idx.length = -1;
   //crop the windows to get the neighbors id
   int miny = MAX(y-radius+1, 1);
   int minx = MAX(x-radius+1, 1);
   int maxy = MIN(y+radius-1, imh);
   int maxx = MIN(x+radius-1, imw);

   int pre_miny = MAX(y-pre_radius+1, 0);
   int pre_minx = MAX(x-pre_radius+1, 0);
   int pre_maxy = MIN(y+pre_radius, imh-1);
   int pre_maxx = MIN(x+pre_radius, imw-1);   
   
   for(i= miny; i< maxy; i++)
       for(j= minx; j< maxx; j++)
       {
           
            if(i>=pre_miny && i<=pre_maxy && j>=pre_minx && j<=pre_maxx)
                continue;
            for(k =0; k<pDim; k++)
            if(img[px3(i, j, k, imh, imw, pDim)]>0){
                
               idx.length = idx.length + 1;
               idx.a[idx.length] = img[px3(i, j, k, imh, imw, pDim)] - 1;
    
           }
       }
   idx.length = idx.length + 1;
   return idx;
}

/*get the cooridinates in the region*/
int* getneighs(int y, int x, double *img, int imh, int imw, int radius, int *numpoint)
{ 
    
   int i,j,k;
   int idx = 0;
   *numpoint = 0;
   
   //crop the windows to get the neighbors id
   int miny = MAX(y-radius+1, 0);
   int minx = MAX(x-radius+1, 0);
   int maxy = MIN(y+radius, imh);
   int maxx = MIN(x+radius, imw);  
   for(i= miny; i< maxy; i++)
    for(j= minx; j< maxx; j++)
    {
        
         *numpoint = *numpoint + img[px2(i,j, imh,ihw)];
   }

   int dim = 2;
   int* subdataPoints= NULL;
   subdataPoints = (int*)malloc((*numpoint)*dim*sizeof(int)); 
   for(i= miny; i<= maxy; i++)
    for(j= minx; j<= maxx; j++)
    {
        for(k=0; k< img[px2(i,j, imh,ihw)]; k++)
        {
            subdataPoints[px2(idx,0,*numpoint,dim)] = i;
            subdataPoints[px2(idx,1,*numpoint,dim)] = j;
            idx = idx + 1;
            
        }
       
   }
   return subdataPoints;
}



