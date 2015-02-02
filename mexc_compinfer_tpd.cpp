#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"        /* the algorithm is connect to matlab */
#include "math.h"
#include "matrix.h"
#include "util.h"
#include "complearn.h"
#include "comptopdown.h"

/* function we use here
 * ===============================================================================
 * Compile: mex -g -I./ mexc_compinfer_tpd.cpp comptopdown.cpp util.cpp filter.cpp checkparts.cpp
 * OR: mex -v CFLAGS="\$CFLAGS -Wall -std=c++0x" -g -I./ mexc_compinfer_tpd.cpp checkparts.cpp complearn.cpp comptopdown.cpp util.cpp filter.cpp
 * 
 * parallel 
 * mex -v CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp -lmpi" COPTIMFLAGS='-O3 -DNDEBUG' -I./ mexc_compinfer_tpd.cpp comptopdown.cpp util.cpp filter.cpp checkparts.cpp
 * 
 * Runing Example: Output = mexc_compinfer_tpd(data,layers, 2, 1,6,10, 0.9, [2 1 1 1 0.5 0.2 0.1], 0.9);
 * ===============================================================================
 *
 * [infer] = mexc_compinfer_tpd(data, layers, numsub, numlayer, numOrient, half_radius, percent, param, discount, rotation);
 *
 */
int totLayers = 7;

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    
    //mtrace();
    
    const mxArray *data;
    const mxArray *layers;
    const mxArray *pinstances;
    const mxArray *ptypes, *pcell, *links, *mask, *morients;
    int numsub;
    int numlayer;
    int numdata;
    int numOrient, rotation;
    int radius;
    double percent;
    double *param;
    double discount;
    int nummax, numel;
    int i,j,k,ifield,subfield, nfields;
    mwIndex jstruct, kstruct;
    mwSize NUM_FIELD_DATA, NUM_OF_DATA, NUM_OF_LAYER, NUM_FIELD_LAYERS;
    mwSize ndim;
    const mwSize *dims;
    mwSize dimsOutput[2];
    double *dyx;
    const mxArray *pprob;
    /*double first, consider to change it into interger*/
    double *ntypes;
    /*definition*/
    const mxArray *elements;
    const mxArray *lambda;
    int numofpi, numoflambda, numofptypes;
    mxClassID datatype;
    // double *elements;
    
    //PROGRAM
    if((nlhs!=1)||(nrhs!=10))
    {
        mexErrMsgTxt("1 output and 10 input needed!\n");
        return;
    }
    
    /*
     * input variable 2: number of subcomposition
    */
    numsub = (int)mxGetScalar(prhs[2]);
    
    /*number of layers */
    numlayer = (int)mxGetScalar(prhs[3]);
    if(numlayer<=0)
    {
        mexErrMsgTxt("Error: The number of Layer should be larger than 0 and less than 7!\n");
        return;
    }
    
    
    /*number of orientations*/
    numOrient  = (int)mxGetScalar(prhs[4]);
    
    /*size of radius*/
    radius  = (int)mxGetScalar(prhs[5]);
    
    /*===========================================
	 * input variable 0: data
     *===========================================
	 */
    data = prhs[0];
    NUM_FIELD_DATA = mxGetNumberOfFields(prhs[0]);
    NUM_OF_DATA =  mxGetNumberOfElements(prhs[0]);
    
    if(NUM_OF_DATA==0)
    {
        mexErrMsgTxt("Error: The data is empty or No input data available!\n");
        return;
    }
    
    dims = mxGetDimensions(data);
    numdata = dims[0] * dims[1];
    
    /*allocation space for the data*/
    struct Data *examples = (Data*)malloc(sizeof(struct Data)*numdata);
    for(i=0; i< numdata; i++) 
    {
        /*allocate numlayer+1 for next layer learning*/
        examples[i].pinstances = (struct sPinstances*)malloc((numlayer)*sizeof(struct sPinstances));
        examples[i].lambda = (struct sLambda *)malloc((numlayer)*sizeof(struct sLambda));
        
    }
    
    
    // M2Map = (const double**)mxCalloc( numdata, sizeof(*M2Map) );   /* MAX1 maps */
    // for(ifield = 0; ifield<NUM_FIELD_LAYERS; ifield++)
    for (jstruct=0; jstruct<NUM_OF_DATA; ++jstruct)
    {


         //int imh = (int)mxGetScalar(mxGetField(prhs[0], jstruct, "imh"));
         //int imw = (int)mxGetScalar(mxGetField(prhs[0], jstruct, "imw"));
         ifield = mxGetFieldNumber(prhs[0],"mask");
         mask =  (const mxArray*)mxGetFieldByNumber(prhs[0], jstruct, ifield);
         //const mxArray *mask =  ( const mxArray *)mxGetPr(mxGetField(prhs[0], jstruct, "mask"));
         dims = mxGetDimensions(mask); 
         int imh = dims[0];
         int imw = dims[1];
         examples[jstruct].mask = (double*)mxGetData(mask);

         //const mxArray *morients =  (const mxArray*)mxGetPr(mxGetField(prhs[0], jstruct, "morients"));
         ifield = mxGetFieldNumber(prhs[0],"morients");
         morients =  ( const mxArray *)mxGetFieldByNumber(prhs[0], jstruct, ifield);
         dims = mxGetDimensions(morients); 
         imh = dims[0];
         imw = dims[1];
         datatype = mxGetClassID(morients);
         if (datatype != mxDOUBLE_CLASS)
         {
             mexErrMsgTxt("warning !! double precision required.");
         }
         examples[jstruct].morients = (double*)mxGetPr(morients);//(unsigned char*)mxGetPr(morients);//
         examples[jstruct].imh = imh;
         examples[jstruct].imw = imw;
         
         /*
           for(i=0; i<dims[0]*dims[1]; i++)
           {

               model[jstruct].ptypes_[numsub-1].p.a[i] = pprob[i];
           }
           model[jstruct].ptypes_[numsub-1].p.length = dims[0]*dims[1];
        */

         ifield = mxGetFieldNumber(prhs[0],"imname");
         if(ifield>-1)
         {
             char *imname =  ( char *)mxArrayToString(mxGetFieldByNumber(prhs[0], jstruct, ifield));
            //char *imname =  (char *)mxGetPr(mxGetField(prhs[0], jstruct, "imname"));
            examples[jstruct].imname = imname;
         }

        examples[jstruct].numlayer = 0;

        /*tmp = mxGetFieldByNumber(prhs[0], jstruct, ifield); */
        // pinstances = (const mxArray *)mxGetPr(mxGetField(prhs[0], jstruct, "pinstances"));
        ifield = mxGetFieldNumber(prhs[0],"pinstances");
        if(ifield<0)
        {    
            examples[jstruct].pinstances = NULL;
            continue;
        }
        pinstances = mxGetFieldByNumber(prhs[0], jstruct, ifield);
        if(pinstances == NULL) {
            mexPrintf("%s%d\t%s%d\n", "FIELD: ", ifield+1, "STRUCT INDEX :", jstruct+1);
            mexErrMsgIdAndTxt( "MATLAB:phonebook:fieldEmpty",
                    "Above field is empty!");
        }
        else
        {
           subfield = mxGetFieldNumber(pinstances, "elements");
           if(subfield<0)
           {    
               examples[jstruct].pinstances[0].elements=NULL;
               continue;
           }
           numofpi = MIN(numlayer, mxGetNumberOfElements(pinstances));
           
           examples[jstruct].numlayer = numofpi;
           //mexPrintf("dataid: %d, pinstance id: %d\n", jstruct, numofpi);
           for(kstruct=0; kstruct <numofpi; kstruct++)
           {
               /*read parts cell*/ 
               // elements = mxGetPr(mxGetField(pinstances, kstruct, "elements"));
               //elements =(const mxArray *)mxGetPr(mxGetField(pinstances, kstruct, "elements"));
               elements = mxGetFieldByNumber(pinstances, kstruct, subfield);
               if(elements==NULL)
               {
                   examples[jstruct].pinstances[kstruct].row = 0;
                   // access the data by point to it
                   examples[jstruct].pinstances[kstruct].col = 0;
                   examples[jstruct].pinstances[kstruct].elements = NULL; 
                   continue;
               }
               else{
                   dims = mxGetDimensions(elements);
                   numel = dims[0] * dims[1];
                  
                   examples[jstruct].pinstances[kstruct].row = dims[0];
                   // access the data by point to it
                   examples[jstruct].pinstances[kstruct].col = dims[1];
                   examples[jstruct].pinstances[kstruct].elements = (double*)mxGetData(elements); 
      
                }
           }
        }


        /*read lambda */

        /*tmp = mxGetFieldByNumber(prhs[0], jstruct, ifield); */
        ifield = mxGetFieldNumber(prhs[0],"lambda");
        if(ifield<0)
            continue;
        //lambda = (const mxArray *)mxGetPr(mxGetField(prhs[0], jstruct, "lambda"));
        lambda = mxGetFieldByNumber(prhs[0], jstruct, ifield);
        if(lambda == NULL) {
            mexPrintf("%s%d\t%s%d\n", "FIELD: ", ifield+1, "STRUCT INDEX :", jstruct+1);
            mexErrMsgIdAndTxt( "MATLAB:phonebook:fieldEmpty",
                    "Above field is empty!");
        }
        else
        {
           /*read parts cell*/ 
           subfield = mxGetFieldNumber(lambda, "elements");
           if(subfield<0)
                continue;
           numoflambda = MIN(numlayer, mxGetNumberOfElements(lambda)); 
           //examples[jstruct].numlayer = numoflambda;
           for(kstruct=0; kstruct <numoflambda; kstruct++)
           {
               // pcell = (const mxArray *)mxGetPr(mxGetField(lambda, kstruct, "elements"));
               pcell = (const mxArray *)mxGetFieldByNumber(lambda, kstruct, subfield);
               if(pcell==NULL){
                   examples[jstruct].lambda[kstruct].elements=NULL;
                   continue;
               }
               dims = mxGetDimensions(pcell);
               int ncells = dims[0] * dims[1];
               /*allocate space for it*/
               examples[jstruct].lambda[kstruct].elements = (cellData*)malloc(ncells*sizeof(cellData));
               for (i=0; i<ncells; ++i)
               {
                    const mxArray* f = mxGetCell(pcell, i);
                    datatype = mxGetClassID(f);
                    if (datatype != mxDOUBLE_CLASS)
                        mexErrMsgTxt("warning !! single precision required.");

                    dyx = (double*)mxGetPr(f);    /* get the pointer to cell content */
                    int h = mxGetM(f);    /* overwriting is ok, since it is constant */
                    int w = mxGetN(f);
                    /*for(j=0;j < w*h; j++)
                        examples[jstruct].lambda[kstruct].elements[i].a[j] = dyx[j];
                     */
                    if(dyx!=NULL){
                        examples[jstruct].lambda[kstruct].elements[i].a = dyx;
                        /*get the length for dyx 8 , 13 , etc*/
                        examples[jstruct].lambda[kstruct].elements[i].length = w*h;
                    }
                    else
                        examples[jstruct].lambda[kstruct].elements[i].a = NULL;


               }
               examples[jstruct].lambda[kstruct].numel = ncells;
           }
        }


    }
    
    
    /*========================================
     * input variable 1: layers / model here
     *========================================
     */
    Parts* model = (Parts*)malloc(totLayers*sizeof(struct Parts));
    
    layers = prhs[1];
    NUM_FIELD_LAYERS = mxGetNumberOfFields(prhs[1]);
    NUM_OF_LAYER =  mxGetNumberOfElements(prhs[1]);
    if(numlayer>NUM_OF_LAYER-1)
    {
        mexPrintf("the current layer %d, is larger than the total number:%d of layer.\n", numlayer, NUM_OF_LAYER);
    }
    /*allocate numlayer+1 for next layer learning*/
    for(jstruct=0; jstruct<NUM_OF_LAYER; jstruct++)
    {
        model[jstruct].ptypes_ = (pTypes*)malloc(sizeof(pTypes)*numsub);
        model[jstruct].ntypes_ = (int*)malloc(sizeof(int)*numsub);
    }
   
    //for(ifield = 0; ifield<NUM_FIELD_LAYERS; ifield++)
        for(jstruct=0; jstruct<NUM_OF_LAYER; jstruct++)
        {
        
            /*tmp = mxGetFieldByNumber(prhs[0], jstruct, ifield); */
        
            ifield = mxGetFieldNumber(prhs[1],"ptypes");
            ptypes =  ( const mxArray *)mxGetFieldByNumber(prhs[1], jstruct, ifield);
            //ptypes = (const mxArray *)mxGetPr(mxGetField(prhs[1], jstruct, "ptypes"));
            if(ptypes == NULL && numlayer !=1) {
                mexPrintf("%s%d\t%s%d\n", "FIELD: ", ifield+1, "STRUCT INDEX :", jstruct+1);
                mexErrMsgIdAndTxt( "MATLAB:phonebook:fieldEmpty",
                        "Above field is empty!");
            }
            else if(ptypes!=NULL)
            {
                
               numofptypes = MIN(numsub, mxGetNumberOfElements(ptypes));  
               for(kstruct=0; kstruct <numofptypes; kstruct++)
               { 
                   /*read parts cell*/ 
                   ifield = mxGetFieldNumber(ptypes, "l");
                   pcell =  ( const mxArray *)mxGetFieldByNumber(ptypes, kstruct, ifield);
                   // pcell = (const mxArray *)mxGetPr(mxGetField(ptypes, numsub-1, "l"));
                   if(pcell==NULL)
                       continue;
                   dims = mxGetDimensions(pcell);
                   int ncells = dims[0] * dims[1];
                   for (i=0; i<ncells; ++i)
                   {
                        const mxArray *f = mxGetCell(pcell, i);
                        datatype = mxGetClassID(f);
                        if (datatype != mxDOUBLE_CLASS)
                            mexErrMsgTxt("warning !! single precision required.");
                        dyx = (double*)mxGetPr(f);    /* get the pointer to cell content */
                        int h = mxGetM(f);    /* overwriting is ok, since it is constant */
                        int w = mxGetN(f);
                        for(int j=0;j < w*h; j++)
                            model[jstruct].ptypes_[kstruct].dyx[i].a[j] = dyx[j];
                        /*get the length for dyx 8 , 13 , etc*/
                        model[jstruct].ptypes_[kstruct].dyx[i].length = w*h;
                   }
                   /*remember the number of parts in each subcomposition*/
                   model[jstruct].ptypes_[kstruct].ntypes = ncells;
                   /*read part probability, array*/
                   ifield = mxGetFieldNumber(ptypes, "p");
                   pprob =  ( const mxArray *)mxGetFieldByNumber(ptypes, kstruct, ifield);
                   //pprob = (const mxArray *)mxGetPr(mxGetField(ptypes, numsub-1, "p"));
                   dims = mxGetDimensions(pprob); 
                   double* temp = (double*)mxGetPr(pprob);
                   for(i=0; i<dims[0]*dims[1]; i++)
                   {
                       model[jstruct].ptypes_[kstruct].p.a[i] = temp[i]; 
                   }
                   model[jstruct].ptypes_[kstruct].p.length = dims[0]*dims[1];
               }
            }
            else
            {
                model[jstruct].ptypes_[kstruct].ntypes = 0;
                printf("Empty part types at the layer: %d\n", jstruct);
            };
            /*get the number of parts*/
            ifield = mxGetFieldNumber(prhs[1],"ntypes");
            ntypes =  ( double*)mxGetPr(mxGetFieldByNumber(prhs[1], jstruct, ifield));
            if(ntypes !=NULL)
            {
                dims = mxGetDimensions(mxGetFieldByNumber(prhs[1], jstruct, ifield));
            
                for(i=0;i<dims[1]*dims[0];i++)
                    model[jstruct].ntypes_[i] = (int)ntypes[i];
                // ntypes = (double *)mxGetPr(mxGetField(prhs[1], jstruct, "ntypes"));
            }
            ifield = mxGetFieldNumber(prhs[1],"links");
            links =  ( const mxArray *)mxGetFieldByNumber(prhs[1], jstruct, ifield);  
            
            //links = (const mxArray *)mxGetPr(mxGetField(prhs[1], jstruct, "links"));
            if(links == NULL) {
                /*mexPrintf("%s%d\t%s%d\n", "FIELD: ", ifield+1, "STRUCT INDEX :", jstruct+1);
                mexErrMsgIdAndTxt( "MATLAB:phonebook:fieldEmpty",
                        "Above field is empty!");
                 */
                 printf("Empty links at the layer: %d\n", jstruct);
                 model[jstruct].plinks_ = NULL;
            }
            else
            {
               /*read parts cell*/ 
               // pcell = (const mxArray *)mxGetPr(mxGetField(links, numsub-1, "l"));
               dims = mxGetDimensions(links);
               int ncells = dims[0] * dims[1];
               //create links 
               model[jstruct].plinks_ = (struct array1d *)malloc(ncells*sizeof(array1d));
               // new array1d[ncells];
               for (i=0; i<ncells; ++i)
               {
                    const mxArray *f = mxGetCell(links, i);
                    datatype = mxGetClassID(f);
                    if (datatype != mxDOUBLE_CLASS)
                        mexErrMsgTxt("warning !! single precision required.");
                    double *bupids = (double*)mxGetPr(f);    /* get the pointer to cell content */
                    int h = mxGetM(f);    /* overwriting is ok, since it is constant */
                    int w = mxGetN(f);
                    //int *temp = bupids;
                    for(j=0;j < w*h; j++)
                        model[jstruct].plinks_[i].a[j] = bupids[j];//links[j];
                    model[jstruct].plinks_[i].length = w*h;
                    
               }
               
            }
            
            ifield = mxGetFieldNumber(prhs[1],"half_radius");
            if(mxGetFieldByNumber(prhs[1], jstruct, ifield))
            model[jstruct].radius =  (int)mxGetScalar(mxGetFieldByNumber(prhs[1], jstruct, ifield));
            // radius = mxGetScalar(mxGetField(prhs[1], jstruct, "half_radius"));
            else model[jstruct].radius = radius;
            
        }//finish read model
        
   
    /*
     * input variable 7, percent to keep the number of parts
     */
    percent = (double)mxGetScalar(prhs[6]);
    
    
    /*
     *input variable 8, average frequency for each layer
     */
    dims = (const int*)mxGetDimensions(prhs[7]);
    ndim = mxGetNumberOfDimensions(prhs[7]);
    mxClassID category;
    category = mxGetClassID(prhs[7]);
    //double* param = mxCreateNumericArray(ndim, dims, category, mxREAL);   
    param=(double *)mxGetData(prhs[7]);       
            
    /*
     * input variable 9: matching discount 
     */
    discount = (double)mxGetScalar(prhs[8]);
     
    /*
     * input variable 10: rotation
     */
    rotation = (int)mxGetScalar(prhs[9]);
    if(rotation >  numOrient)  
    {
            mexPrintf("The input rotation = %d is larger than numOrient = %d\n", rotation, numOrient);
            //mexErrMsgTxt("Error: The number of Layer should be larger than 0 and less than 7!\n");
            return;
    }
    // only 10 parameters in this one

    /*==============================================
     * Inference 
     *==============================================
     */
    struct Data *output;
    //mexPrintf("start inference!\n");
    //output =comptopdown_ex(examples, numdata, model, numlayer, numsub, numOrient,discount, param);
    output =comptopdown_ex2(examples, numdata, model, numlayer, numsub, numOrient,discount, param, rotation);
    //mexPrintf("end of inference!\n");
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    /*
     * output variable 0: inference results
     */
    
    if(NUM_FIELD_LAYERS!=5)
        NUM_FIELD_LAYERS = 5;
    //char *fieldnames[NUM_FIELD_LAYERS]; //This will hold field names.
    char *fieldnames[NUM_FIELD_LAYERS]; //This will hold field names.
    //Assign field names
    for(i=0;i<NUM_FIELD_LAYERS;i++)
        fieldnames[i] = (char*)mxMalloc(20);

    char *subfieldnames[numsub]; //This will hold field names.
    //Assign subfield names
    for(i=0;i<numsub;i++)
        subfieldnames[i] = (char*)mxMalloc(20);
    memcpy(subfieldnames[0],"l",sizeof("l"));
    memcpy(subfieldnames[1],"p",sizeof("p"));
    /*address and size to copy*/
    void *start_of_pr;
    int bytes_to_copy;

    
    memcpy(fieldnames[0],"pinstances",sizeof("pinstances"));
    memcpy(fieldnames[1],"lambda", sizeof("lambda"));
    //Allocate memory for the structure
    plhs[0] = mxCreateStructMatrix(numdata, 1, 2, (const char**)fieldnames);

    memcpy(subfieldnames[0],"elements",sizeof("elements"));
    for(jstruct=0; jstruct<numdata; jstruct++)
    {
        mxArray *elestruct1 = mxCreateStructMatrix(1, 1, 1, (const char**)subfieldnames);
        for(j=0;j<1;j++)
        {   
            ndim = 2;
            /*set dimension for part types: cell*/
            dimsOutput[0]= output[jstruct].pinstances[j].row;
            dimsOutput[1] =output[jstruct].pinstances[j].col;
            /*set dimension for part probability: array*/
            mxArray *instances = mxCreateNumericArray(ndim, dimsOutput, mxDOUBLE_CLASS, mxREAL);
            /* copy
            start_of_pr = (double*)mxGetData(instances);
            bytes_to_copy = dimsOutput[0]*dimsOutput[1]*mxGetElementSize(instances);
            memcpy(start_of_pr, output[jstruct].pinstances[j].elements, bytes_to_copy );
            */
            double *pr= (double*)mxGetData(instances);
            for(int irow=0; irow<dimsOutput[0]; irow++)
                for(int icol=0; icol<dimsOutput[1]; icol++)
                    pr[px2(irow, icol, dimsOutput[0], dimsOutput[1])] = output[jstruct].pinstances[j].elements[px2(irow,icol, dimsOutput[0], dimsOutput[1])];
            //Assign the field values
            mxSetFieldByNumber(elestruct1,j,0, instances);
            //NOTE: mxSetFieldByNumber(..) will automatically take care
        }    
        mxSetFieldByNumber(plhs[0],jstruct,0, elestruct1);
        
        
        
        mxArray *elestruct2 = mxCreateStructMatrix(1, 1, 1, (const char**)subfieldnames);
        for(j=0;j<1;j++)
        {   
            ndim = 2;
            /*set dimension for part types: cell*/
            dimsOutput[0]= output[jstruct].lambda[j].numel;
            dimsOutput[1] = 1;
            mxArray *elemcell = mxCreateCellArray(ndim,dimsOutput);
            
            for(k=0; k< output[jstruct].lambda[j].numel; k++) //
            {
                dimsOutput[0] = output[jstruct].lambda[j].elements[k].length;
                dimsOutput[1] = 2;
                mxArray *f = mxCreateNumericArray(2, dimsOutput, mxDOUBLE_CLASS, mxREAL);
                
                double *pr = (double*)mxGetData(f);
                /*bytes_to_copy = dimsOutput[0]*dimsOutput[1]*mxGetElementSize(f);
                memcpy(start_of_pr, output[jstruct].lambda[j].elements[k].a, bytes_to_copy );
                for(int kk=0;kk< dimsOutput[0]*dimsOutput[1]; kk++)
                    *(pr+kk) = (double)output[jstruct].lambda[j].elements[k].a[kk];
                 **/
                for(int irow=0; irow<dimsOutput[0]; irow++)
                    for(int icol=0; icol<dimsOutput[1]; icol++)
                        pr[px2(irow, icol, dimsOutput[0], dimsOutput[1])] = (double)output[jstruct].lambda[j].elements[k].a[px2(icol,irow, dimsOutput[1], dimsOutput[0])];
                mxSetCell(elemcell, k, f);
            }
            
            //Assign the field values
            mxSetFieldByNumber(elestruct2,j,0, elemcell);
            //NOTE: mxSetFieldByNumber(..) will automatically take care
        }
        mxSetFieldByNumber(plhs[0],jstruct,1, elestruct2);
        
        
    }
    
    /*free subfiled name*/
    for(i=0;i<numsub;i++)
        mxFree(subfieldnames[i]);
    
    /*free the field name*/
    for(i=0;i<NUM_FIELD_LAYERS;i++)
        mxFree(fieldnames[i]);
    
    /*release space for the data*/
    for(i=0; i< numdata; i++) 
    {
        //allocate numlayer+1 for next layer learning//
        if(examples[i].pinstances!=NULL)
        {
            free(examples[i].pinstances);
            examples[i].pinstances = NULL;
        }
        
        for(j=0; j< examples[i].numlayer; j++)
        {
            if(examples[i].lambda[j].elements!=NULL)
            {
                free(examples[i].lambda[j].elements);
                examples[i].lambda[j].elements = NULL;
            }
        }
        
        if(examples[i].lambda!=NULL)
        {
            free(examples[i].lambda);
            examples[i].lambda = NULL;
        }
         
        
    }
    if(examples)
    {
        free(examples);
        examples = NULL;
    }
    
    /*release space for the model*/
    for(jstruct=0; jstruct<NUM_OF_LAYER; jstruct++)
    {
       /*
       for(j=0; j<numsub; j++)
       {
 
           free(model[jstruct].ptypes_[j].dyx);
       }
       */
       free(model[jstruct].ptypes_);
       free(model[jstruct].ntypes_);
       if(model[jstruct].plinks_!=NULL&&jstruct<numlayer)
       {
            free(model[jstruct].plinks_);
            model[jstruct].plinks_ = NULL;
       }
 
    }
    free(model);
    
    //free return value: output
    for(jstruct=0; jstruct<numdata; jstruct++)
    {
        for(j=0;j<1;j++)
        {   
            //mexPrintf("data id: %d, address: %x\n", jstruct, output[jstruct].pinstances[j].elements);
            if(output[jstruct].pinstances[j].elements!= NULL){
                //free it after copy
                free(output[jstruct].pinstances[j].elements);  
                output[jstruct].pinstances[j].elements = NULL;
            }
        }
        // free it 
        if(output[jstruct].pinstances!=NULL)
        {free(output[jstruct].pinstances);output[jstruct].pinstances=NULL;}
        for(j=0;j<1;j++)
        {
            for(k=0; k< output[jstruct].lambda[j].numel; k++){
                if(output[jstruct].lambda[j].elements[k].a!=NULL)
                {
                    //free it after copy
                    free(output[jstruct].lambda[j].elements[k].a);
                    output[jstruct].lambda[j].elements[k].a = NULL;
                }
            }
            //free it
            if(output[jstruct].lambda[j].elements!=NULL)
            {free(output[jstruct].lambda[j].elements);output[jstruct].lambda[j].elements=NULL;}
        }
        // free it
        if(output[jstruct].lambda!=NULL)
        {free(output[jstruct].lambda);output[jstruct].lambda=NULL;}
    }
    if(output!=NULL)
    free(output);
    mexPrintf("done with topdown inference!\n");
    //muntrace();
}
