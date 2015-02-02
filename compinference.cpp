#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "compinference.h"
#include "comptopdown.h"
/*define the total number of layers*/
// int totLayers = 7;
// int numsub = 2;
//#define px2(x, y, lengthx, lengthy) ((x) + (y)*(lengthx) - 0) 
//#define px3(x, y, z, lengthx, lengthy, lengthz) (z+ ((x) + ((y))*(lengthx))*lengthz - 0) 
// struct Parts model[totLayers];

Data* compinference(struct Data *examples, int ndata, struct Parts *model,int numlayer, int numsub,  int numOrient, int radius, 
        double percent, double *param, double discount)
{
    
     Data *output = (Data*)malloc(sizeof(Data)*ndata);
     int i,j,k;
     int imidx, imh,imw;
     int ninst, col, iinst, newid;
     int r, c, o1, o2;
     struct array1d idx1, idx2, idx;
     
     double radjust[] = {0 , 0 , 5, 10, 20 , 30 ,40, 60};
     for(imidx=0; imidx< ndata; imidx++) 
     {
        output[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
        output[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));

     }
     printf("layer: %d, the # of sub-composition: %d, # parts: %d\n",numlayer, numsub, model[numlayer].ntypes_[numsub-1]);
     for(i=0;i<numsub; i++)
         printf("subcomposition: %d,  # parts: %d\n", i, model[numlayer].ntypes_[i]);
     //double **pdist = (double**)malloc(ndata*sizeof(double*));
     int **pdist = (int**)malloc(ndata*sizeof(int*));
     /*inference first*/ 
     //#pragma omp parallel for private(imidx)
     for(imidx=0; imidx< ndata; imidx++)
     {

        int* morients = (int*)examples[imidx].morients;
        imh = examples[imidx].imh;
        imw = examples[imidx].imw;
        ninst = examples[imidx].pinstances[numlayer-1].row;
        col = examples[imidx].pinstances[numlayer-1].col;
        double *pinstances = examples[imidx].pinstances[numlayer-1].elements;
        /*define the output and allocate space*/
        int dimsOutput = 4;
        
        //examples[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
        //examples[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));
        output[imidx].pinstances[0].elements = (double*)malloc(2*ninst*dimsOutput*sizeof(double));
        output[imidx].lambda[0].elements  =  (struct cellData*)malloc(2*ninst*sizeof(struct cellData));
        //initialization the output
        output[imidx].pinstances[0].row = 0;
        output[imidx].pinstances[0].col = 0;
        output[imidx].lambda[0].numel =0;
        if(ninst<1)
        {
            pdist[imidx] = NULL;
            continue;
        }
        //compute the distances between points
        //pdist[imidx] = funDist(pinstances, ninst, col);
        pdist[imidx] = (int*)malloc(ninst*ninst*sizeof(int));
        funDist(pinstances, ninst, col, pdist[imidx]);
        /*new pid in the inference stage*/
        newid = 0;
        int pre_radius = 0;
        //#pragma omp parallel for private(iinst), shared(output)
        for(iinst = 0; iinst< ninst; iinst++)
        {


            r = pinstances[px2(iinst, 1, ninst, col)];
            c = pinstances[px2(iinst, 2, ninst, col)];
            o1 = pinstances[px2(iinst, 0, ninst, col)];
            /* // find the neighbors 
            idx1 = searchneighs(r, c, pinstances, ninst, col, radius, numlayer);
            pre_radius = 4;
            if(numlayer>1)
                pre_radius = model[numlayer-1].radius;
            idx2 = searchneighs2(r, c, pinstances, ninst, col, pre_radius, numlayer-1, radjust);
            // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
            idx = setdiff(idx1,idx2);*/
            pre_radius = 4;
            if(numlayer>1)
                pre_radius = model[numlayer-1].radius - radjust[numlayer-1];
            //idx2 = searchneighs2(r, c, pinstances, ninst, col, pre_radius, numlayer-1, radjust);
            //idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer);  
            //idx = setdiff(idx1,idx2);
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

                    y[i] = pinstances[px2(idx.a[i], 1, ninst, col)];//pinstances[idx.a[i]][2];
                    x[i] = pinstances[px2(idx.a[i], 2, ninst, col)];//pinstances[idx.a[i]][3]; 
                    o2list[i] = pinstances[px2(idx.a[i], 0, ninst, col)];//pinstances[idx.a[i]][1];
                    iy[i] = y[i] - r + radius;
                    ix[i] = x[i] -c + radius;
                    dyx[2*i]= iy[i]- radius;
                    dyx[2*i+1] = ix[i]-radius; 
                    /*if(r==125 && c == 79)
                    printf("the dyx is: %d and %d",y[i], x[i]);*/

                }

                int numdyx = idx.length;
                struct double_array1d ctypeids = matchingparts(o1, r,c, dyx, o2list,numdyx, model, numlayer, numsub, discount);

                int clen = 0;

                if(ctypeids.length>0 && (newid*dimsOutput)<MAX_LEN-1-dimsOutput)
                {
                      for(j=0; j< ctypeids.length; j++)
                      {

                          output[imidx].pinstances[0].elements[newid*dimsOutput + 0] = ctypeids.a[0+5*j]; //pid
                          output[imidx].pinstances[0].elements[newid*dimsOutput + 1] = r;
                          output[imidx].pinstances[0].elements[newid*dimsOutput + 2] = c;
                          output[imidx].pinstances[0].elements[newid*dimsOutput + 3] = ctypeids.a[1 + 5*j];//score
                          /*processing childs*/
                          /*allocate space for the current cell*/
                          output[imidx].lambda[0].elements[newid].a = (double*)malloc(MAX_LEN*sizeof(double));

                          if(examples[imidx].lambda[numlayer-1].numel==0)
                          {
                              output[imidx].lambda[0].elements[newid].length = 0;
                              clen = output[imidx].lambda[0].elements[newid].length;
                              //add itself
                              output[imidx].lambda[0].elements[newid].a[2*clen+0] = r;
                              output[imidx].lambda[0].elements[newid].a[2*clen+1] = c;
                              clen = clen +1;
                              //add children
                              for(k=0; k< ctypeids.a[2 + 5*j]; k++) //for each child
                              {
                                  output[imidx].lambda[0].elements[newid].a[2*clen+0] = y[(int)ctypeids.a[3+ k + 5*j]];
                                  output[imidx].lambda[0].elements[newid].a[2*clen+1] = x[(int)ctypeids.a[3+ k + 5*j]];
                                  clen = clen + 1;
                                  
                                  if(2*clen>MAX_LEN-2)
                                      break;

                              }

                              output[imidx].lambda[0].elements[newid].length = clen;
                          }
                          else
                          {
                              // add itself
                              clen = 0;
                              int len = examples[imidx].lambda[numlayer-1].elements[iinst].length;
                              for(k =0; k< len; k = k+2)
                              {
                                  /*y coordinate*/
                                  output[imidx].lambda[0].elements[newid].a[2*clen+0] = 
                                  examples[imidx].lambda[numlayer-1].elements[iinst].a[k];
                                  /*x coordinate*/
                                  output[imidx].lambda[0].elements[newid].a[2*clen+1] = 
                                  examples[imidx].lambda[numlayer-1].elements[iinst].a[k+1];   
                                  clen = clen + 1;


                              }
                              
                              
                              
                              // copy elements from its children
                              for(k=0; k < ctypeids.a[2 + 5*j]; k++)
                              {    
                                  int chidx = idx.a[(int)ctypeids.a[3 + k + 5*j]]; //instance id

                                  int len = examples[imidx].lambda[numlayer-1].elements[chidx].length;

                                  for(int kk =0; kk< len; kk = kk+2)
                                  {
                                      /*y coordinate*/
                                      output[imidx].lambda[0].elements[newid].a[2*clen+0] = 
                                      examples[imidx].lambda[numlayer-1].elements[chidx].a[kk];
                                      /*x coordinate*/
                                      output[imidx].lambda[0].elements[newid].a[2*clen+1] = 
                                      examples[imidx].lambda[numlayer-1].elements[chidx].a[kk+1];   
                                      clen = clen + 1;

                                  }
                                  
                                  if(2*clen>MAX_LEN-2)
                                        break;

                              }
                              /*final length by accumulate its children*/
                              output[imidx].lambda[0].elements[newid].length = clen;
                          } 

                          newid = newid+1;

                      }
                }

              free(y);
              free(x);
              free(iy);
              free(ix);
              free(o2list);
              free(dyx);

            } 
            
            if(newid>2*ninst-1)
                break;
                
        }
        
        output[imidx].pinstances[0].row = newid;
        output[imidx].pinstances[0].col = dimsOutput;
        output[imidx].lambda[0].numel = newid;
        
        

     }
    //free the allocation for distance computation
    for(imidx=0; imidx<ndata; imidx++)
    {
         if(pdist[imidx]!=NULL)
         {
             free(pdist[imidx]);
             pdist[imidx]=NULL;
         }
    }
    free(pdist);
     
     
     return output;
}


// for cell elements
Data* compinference(struct Data *examples, int ndata, struct Parts *model,int numlayer, int numscale, int numsub,  int numOrient, int radius, 
        double percent, double *param, double discount)
{
    
     Data *output = (Data*)malloc(sizeof(Data)*ndata);
     int i,j,k;
     int imidx, imh,imw;
     int ninst, col, iinst, newid;
     int r, c, o1, o2;
     struct array1d idx1, idx2, idx;
     
     double radjust[] = {0 , 0 , 5, 10, 20 , 30 ,40, 60};
     for(imidx=0; imidx< ndata; imidx++) 
     {
        output[imidx].pinstances = (struct sPinstances*)malloc(numscale*sizeof(struct sPinstances));
        output[imidx].lambda = (struct sLambda *)malloc(numscale*sizeof(struct sLambda));

     }
    
     //double **pdist = (double**)malloc(ndata*numscale*sizeof(double*));
     int **pdist = (int**)malloc(ndata*numscale*sizeof(int*));
     /*inference first*/ 
     int sid=0;//define how many scales
     //#pragma omp parallel for private(imidx)
     for(imidx=0; imidx< ndata; imidx++)
     {
        /*define the output and allocate space*/
        int dimsOutput = 4;
        for(sid=0; sid<numscale; sid++)
        {
            int* morients = (int*)examples[imidx].morients;//did not use at this moment
            imh = examples[imidx].imh;
            imw = examples[imidx].imw;
            ninst = examples[imidx].pinstances[(numlayer-1)*numscale + sid].row;
            col = examples[imidx].pinstances[(numlayer-1)*numscale + sid].col;
            double *pinstances = examples[imidx].pinstances[(numlayer-1)*numscale + sid].elements;
                        
            //examples[imidx].pinstances = (struct sPinstances*)malloc(sizeof(struct sPinstances));
            //examples[imidx].lambda = (struct sLambda *)malloc(sizeof(struct sLambda));
            output[imidx].pinstances[sid].elements = (double*)malloc(2*ninst*dimsOutput*sizeof(double));
            output[imidx].lambda[sid].elements  =  (struct cellData*)malloc(2*ninst*sizeof(struct cellData));
            
            //initialization the output
            output[imidx].pinstances[sid].row = 0;
            output[imidx].pinstances[sid].col = 0;
            output[imidx].lambda[sid].numel =0;
            if(ninst<1)
            {
                pdist[imidx*numscale + sid] = NULL;
                continue;
            }

            //compute the distances between points
            //pdist[imidx*numscale + sid] = funDist(pinstances, ninst, col);
            funDist(pinstances, ninst, col, pdist[imidx*numscale + sid]);
            /*new pid in the inference stage*/
            newid = 0;
            int pre_radius = 0;
            //#pragma omp parallel for private(iinst), shared(output)
            for(iinst = 0; iinst< ninst; iinst++)
            {


                r = pinstances[px2(iinst, 1, ninst, col)];
                c = pinstances[px2(iinst, 2, ninst, col)];
                o1 = pinstances[px2(iinst, 0, ninst, col)];
                /* // find the neighbors 
                idx1 = searchneighs(r, c, pinstances, ninst, col, radius, numlayer);
                pre_radius = 4;
                if(numlayer>1)
                    pre_radius = model[numlayer-1].radius;
                idx2 = searchneighs2(r, c, pinstances, ninst, col, pre_radius, numlayer-1, radjust);
                // idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer); 
                idx = setdiff(idx1,idx2);*/
                pre_radius = 4;
                if(numlayer>1)
                    pre_radius = model[numlayer-1].radius - radjust[numlayer-1];
                //idx2 = searchneighs2(r, c, pinstances, ninst, col, pre_radius, numlayer-1, radjust);
                //idx2 = searchneighs(r, c, pinstances, ninst, col, pre_radius, numlayer);  
                //idx = setdiff(idx1,idx2);
                idx = searchneighs(iinst, pdist[imidx*numscale + sid], ninst, radius, pre_radius, numlayer);
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

                        y[i] = pinstances[px2(idx.a[i], 1, ninst, col)];//pinstances[idx.a[i]][2];
                        x[i] = pinstances[px2(idx.a[i], 2, ninst, col)];//pinstances[idx.a[i]][3]; 
                        o2list[i] = pinstances[px2(idx.a[i], 0, ninst, col)];//pinstances[idx.a[i]][1];
                        iy[i] = y[i] - r + radius;
                        ix[i] = x[i] -c + radius;
                        dyx[2*i]= iy[i]- radius;
                        dyx[2*i+1] = ix[i]-radius; 
                        /*if(r==125 && c == 79)
                        printf("the dyx is: %d and %d",y[i], x[i]);*/

                    }

                    int numdyx = idx.length;
                    struct double_array1d ctypeids = matchingparts(o1, r,c, dyx, o2list,numdyx, model, numlayer, numsub, discount);

                    int clen = 0;

                    if(ctypeids.length>0 && (newid*dimsOutput)<MAX_LEN-1-dimsOutput)
                    {
                          for(j=0; j< ctypeids.length; j++)
                          {

                              output[imidx].pinstances[sid].elements[newid*dimsOutput + 0] = ctypeids.a[0+5*j]; //pid
                              output[imidx].pinstances[sid].elements[newid*dimsOutput + 1] = r;
                              output[imidx].pinstances[sid].elements[newid*dimsOutput + 2] = c;
                              output[imidx].pinstances[sid].elements[newid*dimsOutput + 3] = ctypeids.a[1 + 5*j];//score
                              /*processing childs*/
                              /*allocate space for the current cell*/
                              output[imidx].lambda[sid].elements[newid].a = (double*)malloc(MAX_LEN*sizeof(double));

                              if(examples[imidx].lambda[(numlayer-1)*numscale+sid].numel==0)
                              {
                                  output[imidx].lambda[sid].elements[newid].length = 0;
                                  clen = output[imidx].lambda[sid].elements[newid].length;
                                  //add itself
                                  output[imidx].lambda[sid].elements[newid].a[2*clen+0] = r;
                                  output[imidx].lambda[sid].elements[newid].a[2*clen+1] = c;
                                  clen = clen +1;
                                  //add children
                                  for(k=0; k< ctypeids.a[2 + 5*j]; k++) //for each child
                                  {
                                      output[imidx].lambda[sid].elements[newid].a[2*clen+0] = y[(int)ctypeids.a[3+ k + 5*j]];
                                      output[imidx].lambda[sid].elements[newid].a[2*clen+1] = x[(int)ctypeids.a[3+ k + 5*j]];
                                      clen = clen + 1;

                                      if(2*clen>MAX_LEN-2)
                                          break;

                                  }

                                  output[imidx].lambda[sid].elements[newid].length = clen;
                              }
                              else
                              {
                                  // add itself
                                  clen = 0;
                                  int len = examples[imidx].lambda[(numlayer-1)*numscale + sid].elements[iinst].length;
                                  for(k =0; k< len; k = k+2)
                                  {
                                      /*y coordinate*/
                                      output[imidx].lambda[sid].elements[newid].a[2*clen+0] = 
                                      examples[imidx].lambda[(numlayer-1)*numscale+sid].elements[iinst].a[k];
                                      /*x coordinate*/
                                      output[imidx].lambda[sid].elements[newid].a[2*clen+1] = 
                                      examples[imidx].lambda[(numlayer-1)*numscale+sid].elements[iinst].a[k+1];   
                                      clen = clen + 1;


                                  }



                                  // copy elements from its children
                                  for(k=0; k < ctypeids.a[2 + 5*j]; k++)
                                  {    
                                      int chidx = idx.a[(int)ctypeids.a[3 + k + 5*j]]; //instance id

                                      int len = examples[imidx].lambda[(numlayer-1)*numscale+sid].elements[chidx].length;

                                      for(int kk =0; kk< len; kk = kk+2)
                                      {
                                          /*y coordinate*/
                                          output[imidx].lambda[sid].elements[newid].a[2*clen+0] = 
                                          examples[imidx].lambda[(numlayer-1)*numscale+sid].elements[chidx].a[kk];
                                          /*x coordinate*/
                                          output[imidx].lambda[sid].elements[newid].a[2*clen+1] = 
                                          examples[imidx].lambda[(numlayer-1)*numscale+sid].elements[chidx].a[kk+1];   
                                          clen = clen + 1;

                                      }

                                      if(2*clen>MAX_LEN-2)
                                            break;

                                  }
                                  /*final length by accumulate its children*/
                                  output[imidx].lambda[sid].elements[newid].length = clen;
                              } 

                              newid = newid+1;

                          }
                    }

                  free(y);
                  free(x);
                  free(iy);
                  free(ix);
                  free(o2list);
                  free(dyx);

                } 

                if(newid>2*ninst-1)
                    break;

            }
        
        output[imidx].pinstances[sid].row = newid;
        output[imidx].pinstances[sid].col = dimsOutput;
        output[imidx].lambda[sid].numel = newid;
        
        }//end of each scale

     }
    //free the allocation for distance computation
    for(imidx=0; imidx<ndata; imidx++)
    {
         for(sid = 0; sid<numscale; sid++)
         if(pdist[imidx*numscale + sid]!=NULL)
         {
             free(pdist[imidx*numscale + sid]);
             pdist[imidx*numscale + sid]=NULL;
         }
    }
    free(pdist);
      
    return output;
}



struct double_array1d matchingparts(int o1, int r, int c, int *dyx, int *o2list,int numdyx,
        struct Parts * model, int numlayer,int numsub, double discount)
{
    
  int i,j,k,pid;
  int epsilon = 1.2;  
  //[instances, ctypeids, subpids] = 
  struct double_array1d ctypeids;
  int numctype = 0; //remember the number of matched
  double maxprob = 0;// for best match for a given pid
  
  /*find the typeid index that matching o1*/
  // struct array1d reidx = matchmid(o1, model,numlayer, numsub); 
  /*get the links here:  */ 
  struct array1d reidx = model[numlayer-1].plinks_[o1-1];
  int childid;
  int numchild=0;
  int nummatched =0;
  double score = 0;
  
  int pDim = 3;
  pDim = pDim + MIN(numsub,2); //because the maximum number of child is 2
  
  for(i=0; i< reidx.length; i++)
  {
     pid = reidx.a[i]-1;//start from 0
     numchild = model[numlayer].ptypes_[numsub-1].dyx[pid].a[2];
     int *flag = (int*)malloc(sizeof(int)*numchild);
     int *bestk = (int*)malloc(sizeof(int)*numchild); // for the best id in dyx vector
     score = 1;
     /* match two childrens */
     for(j=0; j< numchild; j = j+1)
     {   
         
         flag[j] = 0; //initialize the matching flag
         childid = model[numlayer].ptypes_[numsub-1].dyx[pid].a[3+5*j];
         //double* mean = (double*)malloc(2*sizeof(double));
         //double* var = (double*)malloc(2*sizeof(double));
         double mean[2], var[2];
         /*find the matched child orientation */
         int* indicator = find(o2list, numdyx, childid, 'e');
         /*THE NUM OF MATCHED*/
         nummatched = sum(indicator, numdyx);
         maxprob = 0;// for best match for a given pid
         if(nummatched>0) 
         {
             //double* selected_dyx = selection(dyx, numdyx, indicator);
             //double* selected_o2 = selection(o2list, numdyx, indicator);
             //then match relative positions
             mean[0] = model[numlayer].ptypes_[numsub-1].dyx[pid].a[4+5*j];
             mean[1] = model[numlayer].ptypes_[numsub-1].dyx[pid].a[5+5*j];
             var[0] = model[numlayer].ptypes_[numsub-1].dyx[pid].a[6+5*j];
             var[1] = model[numlayer].ptypes_[numsub-1].dyx[pid].a[7+5*j];
             
             for(k=0; k<numdyx; k++)
             {
                 if(indicator[k]==1)
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

             if(maxprob > 0.2*0.5*pow(discount,numlayer-1))
             { //it means a match
                flag[j] = 1;
                score = score*maxprob;
             }
             /*
             for(k=0; k<numdyx; k++)
             {
                 if(indicator[k]==1)
                 {
                     double prob = (0.5*(dyx[2*k]- mean[0])*(dyx[2*k]- mean[0]))/(var[0]+pow(epsilon, numlayer-1)) + 0.5*((dyx[2*k+1]- mean[1])*(dyx[2*k+1]- mean[1]))/(var[1]+pow(epsilon,numlayer-1)) ;
                     maxprob = exp(-prob);
                     if(maxprob > 0.1*pow(discount,numlayer-1))//if(maxprob > 0.1*pow(discount,numlayer-1))
                     {
                        
                        flag[j] = 1;
                        score = score*maxprob;
                        bestk[j] = k;
                        break;
                     }
                 }
             }*/

         }
         if(indicator)
         {free(indicator);indicator=NULL;}
         //free(mean);free(var);
     }
     
     /*if both of the children match, then it has a matching here*/
     if(sum(flag, numchild) == numchild)
     { //it means a match
        //remember the matched typeid in the array [olist and dyx]  
        ctypeids.a[numctype*pDim+0] = pid+1;//add back to the original pid
        ctypeids.a[numctype*pDim+1] = pow(score, 1/((double)numchild));
        //printf("matched score: %f\n", ctypeids.a[numctype*pDim+1]);
        ctypeids.a[numctype*pDim+2] = numchild;
        for(j=0; j< numchild; j = j+1)
            ctypeids.a[numctype*pDim+3+j] = bestk[j]; //% store the matched id in dyx
        
        numctype = numctype + 1;
        if(numctype*pDim>MAX_LEN-pDim-1)
            break;
     }
     free(flag);
     free(bestk);
     
  } 
  ctypeids.length = numctype;
  if(ctypeids.length<2) //no match or one match
      return ctypeids;
  
  //printf("detected childrens: %d\n", numctype);
  /*how many two-subcomposition*/
  nummatched = 0;
  int* childlen = (int*)malloc(sizeof(int)*numctype);
  /* consider suppression here */
  for(i=0;i<numctype; i++)
  {
      if(ctypeids.a[i*pDim+2]==2) //consider two children first
            nummatched = nummatched + 1;
      //compute how many childs for each pid here
      //childlen[i] = getnumchild(ctypeids.a[i*pDim+0], model, numlayer+1, numsub);
  }
  
  int idx=0;
  if(nummatched==0 && numctype>0)
  {
      score = 0;
      double score2 = 0;
      int idx2 =0;
      for(i=0;i<numctype; i++){
          //if(score<childlen[i])//ctypeids.a[i*pDim+1])
          if(score< ctypeids.a[i*pDim+1])
          {
               score2 = score;
               idx2 = idx;
               score = ctypeids.a[i*pDim+1];
               //score = childlen[i];
               idx = i;
          }else{ 
              if(score2< score && score2 < ctypeids.a[i*pDim+1])
              {
                   score2 = ctypeids.a[i*pDim+1];
                   idx2 = i;

              }
          }
      
      }
      if(numctype==1)
      {
          ctypeids.a[0]= ctypeids.a[idx*pDim+0];//get the best pid
          ctypeids.a[1]= ctypeids.a[idx*pDim+1];//get the best matched score
          ctypeids.a[2]= ctypeids.a[idx*pDim+2];//get the # childs
          for(j=0; j< ctypeids.a[idx*pDim+2];j++)
               ctypeids.a[3+j]= ctypeids.a[idx*pDim+3+j];//get the # childs

          ctypeids.length =1;
      }
      else
      {
          if(idx2==0)
          {
              ctypeids.a[0+pDim]= ctypeids.a[idx*pDim+0];//get the best pid
              ctypeids.a[1+pDim]= ctypeids.a[idx*pDim+1];//get the best matched score
              ctypeids.a[2+pDim]= ctypeids.a[idx*pDim+2];//get the # childs
              for(j=0; j< ctypeids.a[idx*pDim+2];j++)
                   ctypeids.a[3+j+pDim]= ctypeids.a[idx*pDim+3+j];//get the # childs
              ctypeids.length =2;
          }
          else{

              ctypeids.a[0]= ctypeids.a[idx*pDim+0];//get the best pid
              ctypeids.a[1]= ctypeids.a[idx*pDim+1];//get the best matched score
              ctypeids.a[2]= ctypeids.a[idx*pDim+2];//get the # childs
              for(j=0; j< ctypeids.a[idx*pDim+2];j++)
                   ctypeids.a[3+j]= ctypeids.a[idx*pDim+3+j];//get the # childs
              ctypeids.length =2;


              ctypeids.a[0+pDim]= ctypeids.a[idx2*pDim+0];//get the second best pid
              ctypeids.a[1+pDim]= ctypeids.a[idx2*pDim+1];//get the second best matched score
              ctypeids.a[2+pDim]= ctypeids.a[idx2*pDim+2];//get the # childs
              for(j=0; j< ctypeids.a[idx2*pDim+2];j++)
                   ctypeids.a[3+j+pDim]= ctypeids.a[idx2*pDim+3+j];//get the # childs
              ctypeids.length =2;


          }
      }
      
  }
  else if(nummatched>1 && numctype>0) //more than two matches
  {
      score = 0;
      double score2 = 0;
      int idx2 =0;
      for(i=0;i<numctype; i++)
      {
          if(ctypeids.a[i*pDim+2] == 2){
              //if(score<ctypeids.a[i*pDim+1] && ctypeids.a[i*pDim+2] == 2)
              if(score<ctypeids.a[i*pDim+1])
              {
                   score2 = score;
                   idx2 = idx;
                   score = ctypeids.a[i*pDim+1];
                   idx = i;

              }else {
                  if(score2<=ctypeids.a[i*pDim+1])
                  {
                       score2 = ctypeids.a[i*pDim+1];
                       idx2 = i;
                  }
              }
          }
      }
      
      if(idx2==0)
      {
          ctypeids.a[0+pDim]= ctypeids.a[idx*pDim+0];//get the best pid
          ctypeids.a[1+pDim]= ctypeids.a[idx*pDim+1];//get the best matched score
          ctypeids.a[2+pDim]= ctypeids.a[idx*pDim+2];//get the # childs
          for(j=0; j< ctypeids.a[idx*pDim+2];j++)
               ctypeids.a[3+j+pDim]= ctypeids.a[idx*pDim+3+j];//get the # childs
          ctypeids.length =2;
      }
      else{
          
          ctypeids.a[0]= ctypeids.a[idx*pDim+0];//get the best pid
          ctypeids.a[1]= ctypeids.a[idx*pDim+1];//get the best matched score
          ctypeids.a[2]= ctypeids.a[idx*pDim+2];//get the # childs
          for(j=0; j< ctypeids.a[idx*pDim+2];j++)
               ctypeids.a[3+j]= ctypeids.a[idx*pDim+3+j];//get the # childs
          ctypeids.length =2;
          
          
          ctypeids.a[0+pDim]= ctypeids.a[idx2*pDim+0];//get the second best pid
          ctypeids.a[1+pDim]= ctypeids.a[idx2*pDim+1];//get the second best matched score
          ctypeids.a[2+pDim]= ctypeids.a[idx2*pDim+2];//get the # childs
          for(j=0; j< ctypeids.a[idx2*pDim+2];j++)
               ctypeids.a[3+j+pDim]= ctypeids.a[idx2*pDim+3+j];//get the # childs
          ctypeids.length =2;
          
          
      }
  }
  
  else if(nummatched>0 && numctype>0) //one with two children matched
  {
      score = 0;
      idx = 0;
      for(i=0;i< numctype; i++)
          if(score<=ctypeids.a[i*pDim+1] && ctypeids.a[i*pDim+2] == 2)
          //if(score<childlen[i])
          {
               score = ctypeids.a[i*pDim+1];
               idx = i;
          }
      ctypeids.a[0]= ctypeids.a[idx*pDim+0];//get the best pid
      ctypeids.a[1]= ctypeids.a[idx*pDim+1];//get the best matched score
      ctypeids.a[2]= ctypeids.a[idx*pDim+2];//get the # childs
      for(j=0; j< ctypeids.a[idx*pDim+2];j++)
           ctypeids.a[3+j]= ctypeids.a[idx*pDim+3+j];//get the # childs
      ctypeids.length =1;
      
  }
  else
      ctypeids.length = 0;
  
  free(childlen);
  return ctypeids;
    
}


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
            /*save pid here*/
            reidx.a[count] = model[numlayer].ptypes_[numsub-1].dyx[i].a[0];
            count = count + 1;
        }
    
    }
    reidx.length = count;
    return reidx;
    
}



