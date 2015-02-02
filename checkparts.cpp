#include <stdio.h>
#include <math.h>
#include "checkparts.h"
/*matching parts with orientation and relative positions*/
struct array1d checksubparts(int o1, int r, int c, int* dyx, int* o2list, int numdyx, struct pTypes *ptypes, int* nptypes, 
        int numlayer, int subindex, double discount, struct array1d *lutlinks)
{
  int i,j,pid;
  int epsilon = 1.1;  
  //[instances, ctypeids, subpids] = 
  struct array1d ctypeids, sids;
  struct double_array1d pscores;
  int numctype = 0;
  double maxprob = 0;// for best match for a given pid     
  double  mean[2];
  double  var[2];
  /*find the typeid index that matching o1*/
  struct array1d reidx;
  if(lutlinks)
      reidx = lutlinks[o1-1];
  else
      reidx= searchanchorid(o1, ptypes, nptypes, subindex);  
  
  for(i=0; i< reidx.length; i++)
  {
     pid = reidx.a[i];
     //double* dyx = (double*)malloc((pDim-1)*sizeof(double));
     int childid = ptypes[subindex-1].dyx[pid].a[3];
     /*find the matched child orientation */
     int* indicator = find(o2list, numdyx, childid, 'e');
     /*THE NUM OF MATCHED*/
     int nummatched = sum(indicator, numdyx);
     maxprob = 0;// for best match for a given pid
     if(nummatched > 0) 
     {
         //double* selected_dyx = selection(dyx, numdyx, indicator);
         //double* selected_o2 = selection(o2list, numdyx, indicator);
         //then match relative positions
         mean[0] = ptypes[subindex-1].dyx[pid].a[4];//y 
         mean[1] = ptypes[subindex-1].dyx[pid].a[5];//x
         var[0] = ptypes[subindex-1].dyx[pid].a[6];
         var[1] = ptypes[subindex-1].dyx[pid].a[7];
         for(j=0; j<numdyx; j++)
         {
             if(indicator[j]==1)
             {
                 double prob = (0.5*(dyx[2*j]- mean[0])*(dyx[2*j]- mean[0]))/(var[0]+pow(epsilon, numlayer-1)) + 0.5*((dyx[2*j+1]- mean[1])*(dyx[2*j+1]- mean[1]))/(var[1]+pow(epsilon,numlayer-1)) ;
                 //prob = exp(-prob);
                 if(maxprob<exp(-prob))
                     maxprob = exp(-prob);

             }
         }
         
         if(maxprob > 0.5*pow(discount,numlayer-1))
         { 
            
            if(numctype>MAX_LEN-2)
                break;
            //it means a match
            ctypeids.a[numctype] = pid;//remember the matched typeid in the array [olist and dyx]
            pscores.a[numctype] = maxprob;//remember the scores
            numctype = numctype + 1;

         }
     }
     if(indicator)
     {free(indicator);indicator = NULL;}

  } 

  
  ctypeids.length = numctype;
  
  if(ctypeids.length>1)
  {
      double angi, angj;   
      //only consider the parts without angle overlap in numOrient
      int *indicator= (int*)malloc(ctypeids.length*sizeof(int));
      //printf("detected: %d, idlist=", ctypeids.length);
      for(i=0; i< ctypeids.length; i++)
      {      
          indicator[i] = 1;
          //printf("%2.2f ", pscores.a[i]);
      }
      //printf("\n");
      for(i=0; i< ctypeids.length-1; i++)
      {
        for(j=i+1; j< ctypeids.length; j++)
        {
            //first angle
            pid = ctypeids.a[i];
            mean[0] = ptypes[subindex-1].dyx[pid].a[4];//y 
            mean[1] = ptypes[subindex-1].dyx[pid].a[5];//x
            angi = atan2(mean[0],mean[1]);//*180/ PI;
            //angi = atan(mean[0]/mean[1]);
            //printf("pid =%d: y1= %f x1=%f angle = %f, ",pid, mean[0], mean[1], angi);
            //second angle
            pid = ctypeids.a[j];
            mean[0] = ptypes[subindex-1].dyx[pid].a[4];//y 
            mean[1] = ptypes[subindex-1].dyx[pid].a[5];//x
            angj = atan2(mean[0],mean[1]);// * 180 / PI;
            //printf("pid =%d: y2= %f x2=%f angle = %f\n",pid, mean[0], mean[1],angj);
            //check the orientation overlap & suppression
            if(fabs(angj-angi)< double(PI)/4)//0.3927
            {
                //printf("length=%d, pid=%d scorei= %f,pjd=%d scorej= %f, diff= %f threshold=%f\n",ctypeids.length, ctypeids.a[i], pscores.a[i], ctypeids.a[j], pscores.a[j], fabs(angj-angi), double(PI)/8);
                //printf("pid1=%d pid2= %d, diff= %f threshold=%f\n", ctypeids.a[i], ctypeids.a[j], fabs(angj-angi), double(PI)/8);
                if((pscores.a[i]>pscores.a[j])&&indicator[j]>0)
                    indicator[j] = 0;
                else if((pscores.a[i]<pscores.a[j])&&indicator[i]>0)
                    indicator[i] = 0;
                else
                {
                    
                }
            }
        }
      }
      /*
      for(i=0; i< ctypeids.length; i++)
      {      
          printf(" %d ", indicator[i]);
      }
      printf("\n");
      
      if(ctypeids.length==3)
      {
          for(i=0; i< ctypeids.length; i++)
          {      
              
              printf("pid=%d kept=%d ", ctypeids.a[i], indicator[i]);
          }
          printf("\n");
      }
      */
      j = 0;sids.length=0;
      for(i=0; i< ctypeids.length; i++)
      {
          
         if(indicator[i]==1){
            sids.a[j] = ctypeids.a[i];
            j = j+1;
         }
      }
      sids.length = j;
      //printf("original length= %d, after = %d\n",ctypeids.length, sids.length);
      free(indicator);
      return sids;
  }   
  else
  return ctypeids;
    
}

/*find the parts whoes anchorid equals to o1*/
struct array1d searchanchorid(int o1, struct pTypes *ptypes, int *nptypes, int subindex)
{
    struct array1d reidx;
    int count = 0;
    for(int i=0; i< nptypes[subindex-1]; i++)
    {
        
        if(ptypes[subindex-1].dyx[i].a[1] == o1)
        {
        
            reidx.a[count] = i;
            count = count + 1;
   
        }
    
    
    }
    reidx.length = count;
    return reidx;
}

/*suppression straight line*/
double *setweights(struct pTypes *ptypes, int *nptypes, int subindex, double rate)
{
    int i,j;
    int o1, o2;
    double *weights = (double*)malloc(sizeof(double)*nptypes[subindex-1]);
    for(i=0; i< nptypes[subindex-1]; i++)
    {
        weights[i] = 1;
    }
    int count = 0;
    for(i=0; i< nptypes[subindex-1]; i++)
    {   
        o1 = ptypes[subindex-1].dyx[i].a[1];
        for(j=0; j< ptypes[subindex-1].dyx[i].a[2]; j++) 
        {
        
            o2 = ptypes[subindex-1].dyx[i].a[3+ 5*j];
            if(o1 == o2)
                weights[i] = weights[i]*rate;
   
        }
    	if(ptypes[subindex-1].dyx[i].a[2]==2)
        {
            if(ptypes[subindex-1].dyx[i].a[3] == ptypes[subindex-1].dyx[i].a[8])
                weights[i] = weights[i]*rate;
        }
    
    }
    
    return weights;
}

double adjweight(struct pTypes *ptypes, int *nptypes, int subindex, int anchorid, int cid, double rate)
{
    int i,j;
    int o1, o2, o3;
    double temp=1;

    if(anchorid > nptypes[subindex-1]-1 || cid > nptypes[subindex-1]-1)
        printf("index exceeding error: anchor = %d, cid = %d, total number of parts: %d\n", anchorid, cid, nptypes[subindex-1]);
    o1 = ptypes[subindex-1].dyx[anchorid].a[1];
    for(j=0; j< ptypes[subindex-1].dyx[anchorid].a[2]; j++) 
    {

        o2 = ptypes[subindex-1].dyx[anchorid].a[3+ 5*j];
        if(o1 == o2)
            temp = temp*rate;

    }
    
    o3 = ptypes[subindex-1].dyx[cid].a[1];;
    if(o1 == o3 || o1==ptypes[subindex-1].dyx[cid].a[3+ 5*j])
    {   
        for(j=0; j< ptypes[subindex-1].dyx[cid].a[2]; j++) 
        {
            if(ptypes[subindex-1].dyx[cid].a[3+ 5*j] == o2)
                temp = temp*rate;

        }
    }
    else
        printf("anchor matching error: anchor = %d, cid = %d\n", anchorid, cid);
    
    return temp;
    
}



/*suppression straight line*/
double * setweights(double *ptypes, int *nptypes,  int pDim, int subindex, double rate)
{
    int i,j;
    int o1, o2;
    double *weights = (double*)malloc(sizeof(double)*nptypes[subindex-1]);
    for(i=0; i< nptypes[subindex-1]; i++)
    {
        weights[i] = 1;
    }
    int count = 0;
    for(i=0; i< nptypes[subindex-1]; i++)
    {   
        o1 = ptypes[i*pDim+1];
        for(j=0; j< ptypes[i*pDim+2]; j++) 
        {
        
            o2 = ptypes[i*pDim + 3+ 5*j];
            if(o1 == o2)
                weights[i] = weights[i]*rate;
   
        }
    
    
    }
    
    return weights;
}


double adjweight(double *ptypes, int *nptypes, int pDim, int subindex, int anchorid, int cid, double rate)
{
    int i,j;
    int o1, o2, o3;
    double temp=1;

    if(anchorid > nptypes[subindex-1]-1 || cid > nptypes[subindex-1]-1)
        printf("index exceeding error: anchor = %d, cid = %d, total number of parts: %d\n", anchorid, cid, nptypes[subindex-1]);
    o1 = ptypes[anchorid*pDim+1];
    for(j=0; j< ptypes[anchorid*pDim+2]; j++) 
    {

        o2 = ptypes[anchorid*pDim + 3+ 5*j];
        if(o1 == o2)
            temp = temp*rate;

    }
    
    o3 = ptypes[cid*pDim+1];
    if(o1 == o3 || o1==ptypes[cid*pDim + 3+ 5*j])
    {   
        for(j=0; j< ptypes[cid*pDim+2]; j++) 
        {
            if(ptypes[cid*pDim + 3+ 5*j] == o2)
                temp = temp*rate;

        }
    }
    else
        printf("anchor matching error: anchor = %d, cid = %d\n", anchorid, cid);
    
    return temp;
    
}



struct array1d checksubparts(int o1, int r, int c, int* dyx, int* o2list, int numdyx, double* ptypes, 
        int* nptypes, int pDim,int numlayer, int subindex, double discount, struct array1d *lutlinks)
{
  int i,j,pid;
  int epsilon = 2.5;  
  //[instances, ctypeids, subpids] = 
  struct array1d ctypeids;
  int numctype = 0;
  double maxprob = 0;// for best match for a given pid
  /*find the typeid index that matching o1*/
  struct array1d reidx;
  if(!lutlinks)
      reidx = lutlinks[o1-1];
  else
      reidx= searchanchorid(o1, ptypes, nptypes, pDim);  
  for(i=0; i< reidx.length; i++)
  {
     pid = reidx.a[i];
     //double* dyx = (double*)malloc((pDim-1)*sizeof(double));
     int childid = ptypes[pid*pDim+3];
     double* mean = (double*)malloc(2*sizeof(double));
     double* var = (double*)malloc(2*sizeof(double));
     /*find the matched child orientation */
     int* indicator = find(o2list, numdyx, childid, 'e');
     /*THE NUM OF MATCHED*/
     int nummatched = sum(indicator, numdyx);
     maxprob = 0;// for best match for a given pid
     if(nummatched > 0) 
     {
         //double* selected_dyx = selection(dyx, numdyx, indicator);
         //double* selected_o2 = selection(o2list, numdyx, indicator);
         //then match relative positions
         mean[0] = ptypes[pid*pDim+4];
         mean[1] = ptypes[pid*pDim+5];
         var[0] = ptypes[pid*pDim+6];
         var[1] = ptypes[pid*pDim+7];
         for(j=0; j<numdyx; j++)
         {
             if(indicator[j]==1)
             {
                 double prob = (0.5*(dyx[2*j]- mean[0])*(dyx[2*j]- mean[0]))/(var[0]+pow(epsilon, numlayer)) + 0.5*((dyx[2*j+1]- mean[1])*(dyx[2*j+1]- mean[1]))/(var[1]+pow(epsilon,numlayer)) ;
                 //prob = exp(-prob);
                 if(maxprob<exp(-prob))
                     maxprob = exp(-prob);

             }
         }
         
         
         if(maxprob > 0.5*pow(discount,numlayer-1))
         { //it means a match
            ctypeids.a[numctype] = pid;//remember the matched typeid in the array [olist and dyx]
            numctype = numctype + 1;
            if(numctype>MAX_LEN-1)
                break;
         }
         
     }
     free(indicator);
     free(mean);
     free(var);

  } 

  ctypeids.length = numctype;
  return ctypeids;
    
}


struct array1d searchanchorid(int o1, double* ptypes, int *nptypes, int pDim)
{
    struct array1d reidx;
    int count = 0;
    for(int i=0; i< nptypes[0]; i++)
    {
        
        if(ptypes[pDim*i+1] == o1)
        {
        
            reidx.a[count] = i;
            count = count + 1;
   
        }
    
    
    }
    reidx.length = count;
    return reidx;
}



/* selection */
double* selection(int* dyx, int numdyx, int dyxDim, int* indicator)
{
    
    int nomatched = sum(indicator, numdyx);
    double* selected_dyx = (double*)malloc(nomatched*dyxDim*sizeof(double)); 
    
    for(int i=0; i< nomatched; i++)
    {
        
    
    }
    return selected_dyx;
}
