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
 * Compile: mex -g -I./ mexc_comptopdown.cpp comptopdown.cpp util.cpp filter.cpp checkparts.cpp
 * Runing Example: Output = mexc_comptopdown(data,layers, 2, 5, 4, 0.9, [2 1 1 1 0.5 0.2 0.1]);
 * ===============================================================================
 *
 comptoptopdown(int *morients, int imh, int imw, int layers, int numlayer,int numOrient,int numsub, double discount,double *param)
 */
int totLayers = 7;

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    
    //mtrace();
    
    const mxArray *data;
    /*
    const mxArray *layers;
    const mxArray *pinstances;*/
    const mxArray *ptypes, *pcell, *links;// *mask, *morients;
    int numsub;
    int numlayer;
    int numdata;
    int numOrient;
    int *radius;
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
    int numofpi;
    mxClassID datatype;
    // double *elements;
    bool flag = false;
    //PROGRAM
    if((nlhs!=1)||(nrhs!=7))
    {
        mexErrMsgTxt("1 output and 9 input needed!\n");
        return;
    }
    
    
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
    //numdata = dims[0] * dims[1];
    // image size
    int imh = dims[0];
    int imw = dims[1];
    double *morients = (double*)mxGetData(data);

    //set it for testing only
    numdata =1;
    
    /*
     * input variable 2: number of subcomposition
    */
    numsub = (int)mxGetScalar(prhs[2]);
    if(numsub!=2)
    {
        mexErrMsgTxt("Error: only consider 2-sub composition!\n");
        return;
        
    }
    /*number of layers */
    numlayer = (int)mxGetScalar(prhs[3]);
    if(numlayer<=0)
    {
        mexErrMsgTxt("Error: The number of Layer should be larger than 0 and less than 7!\n");
        return;
    }
    
    
    /*number of orientations*/
    numOrient  = (int)mxGetScalar(prhs[4]);  
    
    
    
    /*
     * input variable 6: matching discount 
     */
    discount = (double)mxGetScalar(prhs[5]);
    
        /*
     *input variable 7, parameters for each layer
     */
    param = (double*)mxGetData(prhs[6]);
    dims = (const int*)mxGetDimensions(prhs[6]);
    ndim = mxGetNumberOfDimensions(prhs[6]);
    mxClassID category;
    category = mxGetClassID(prhs[6]);
    //double* param = mxCreateNumericArray(ndim, dims, category, mxREAL);   
    param=(double *)mxGetData(prhs[6]);       
    
    // only 7 parameters in this one
    
    
    /*========================================
     * input variable 1: layers / model here
     *========================================
     */
    Parts* model = (Parts*)malloc(totLayers*sizeof(struct Parts));
    
    //layers = prhs[1];
    NUM_FIELD_LAYERS = mxGetNumberOfFields(prhs[1]);
    NUM_OF_LAYER =  mxGetM(prhs[1]);
    if(numlayer>NUM_OF_LAYER){
        mexErrMsgTxt("the current layer is larger than the total number of layer.\n");
        return;
    }
    /*allocate numlayer+1 for next layer learning*/
    for(jstruct=0; jstruct<numlayer; jstruct++)
    {
        model[jstruct].ptypes_ = (pTypes*)malloc(sizeof(pTypes)*(numsub+1));
        model[jstruct].ntypes_ = (int*)malloc(sizeof(int)*(numsub+1));
        
        for(j=0;j<numsub+1; j++)
        {    
            model[jstruct].ptypes_[j].ntypes=0;
            model[jstruct].ntypes_[j]=0;
        }
    }
   
    for(jstruct=0; jstruct<numlayer; jstruct++)
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

           //numofpi = MIN(numsub, mxGetNumberOfElements(ptypes));
           if(mxGetNumberOfElements(ptypes)>numsub)
                    numofpi = numsub+1;  
               else
                   numofpi = numsub; 
           for(kstruct=0; kstruct <numofpi; kstruct++)
           { 
               /*read parts cell*/ 
               ifield = mxGetFieldNumber(ptypes, "l");
               pcell =  ( const mxArray *)mxGetFieldByNumber(ptypes, kstruct, ifield);
               // pcell = (const mxArray *)mxGetPr(mxGetField(ptypes, numsub-1, "l"));
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
           //consider the input format for ptypes
           
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
        model[jstruct].radius =  (int)mxGetScalar(mxGetFieldByNumber(prhs[1], jstruct, ifield));
        
    }//finish read model
    
    //consider extra format in the code
    if((mxGetNumberOfElements(ptypes)>numsub)&&model[numlayer-1].ptypes_[numsub].ntypes>0)//only the last layers
        flag = true;
    /*if(flag)
        numsub = numsub + 1;*/
    /*==============================================
     * Inference 
     *==============================================
     */
    struct Data *output;
    output =comptopdown(morients,imh, imw, model, numlayer, numsub, numOrient, discount, param);
    
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
    
//     /*release space for the data*/
//     for(i=0; i< numdata; i++) 
//     {
//         /*allocate numlayer+1 for next layer learning*/
//         if(examples[i].pinstances)
//         {
//             free(examples[i].pinstances);
//             examples[i].pinstances = NULL;
//         }
//         
//         for(j=0; j< numlayer; j++)
//         {
//             if(examples[i].lambda[j].elements)
//             {
//                 free(examples[i].lambda[j].elements);
//                 examples[i].lambda[j].elements = NULL;
//             }
//         }
//         
//         if(examples[i].lambda)
//         {
//             free(examples[i].lambda);
//             examples[i].lambda = NULL;
//         }
//          
//         
//     }
//     if(examples)
//     {
//         free(examples);
//         examples = NULL;
//     }
    
    /*release space for the model*/
    for(jstruct=0; jstruct<numlayer; jstruct++)
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
            //free it after copy
            free(output[jstruct].pinstances[j].elements);  
        }
        // free it 
        free(output[jstruct].pinstances);
        for(j=0;j<1;j++)
        {
            for(k=0; k< output[jstruct].lambda[j].numel; k++)
                //free it after copy
                free(output[jstruct].lambda[j].elements[k].a);
            //free it
            free(output[jstruct].lambda[j].elements);
        }
        // free it
        free(output[jstruct].lambda);
    }
    free(output);
    
    mexPrintf("done!\n");
    //muntrace();
}
