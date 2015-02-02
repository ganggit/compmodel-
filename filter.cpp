#include "filter.h"
void createFilter(int radius, double sigma, double ** &gKernel)
{   
    
    // gKernel = double_matrix(2*klen+1, 2*klen+1);
    gKernel = double_matrix(2*radius+1, 2*radius+1);
    // set standard deviation to 1.0
    //double sigma = 1.0;
    double r, s = 2.0 * sigma * sigma;
 
    // sum is for normalization
    double sum = 0.0;
 
    // generate 5x5 kernel
    for (int x = -radius; x <= radius; x++)
    {
        for(int y = -radius; y <= radius; y++)
        {
            r = sqrt(x*x + y*y);
            gKernel[x + radius][y+radius] = (exp(-(r*r)/s))/(PI * s);
            sum += gKernel[x + radius][y + radius];
        }
    }
 
    // normalize the Kernel
    for(int i = 0; i < 2*radius+1; ++i)
        for(int j = 0; j < 2*radius+1; ++j)
            gKernel[i][j] /= sum;
 
}


double* conv2(double* in, int rows, int cols, double **kernel,int kRows, int kCols)
{
    int i,j,m,n;
    int ii,jj,mm,nn;
    double sum =0;
    // find center position of kernel (half of kernel size)
    int kCenterX = round((kCols-1)/2);
    int kCenterY = round((kRows-1)/2);
    double* out = (double*)malloc(sizeof(double)*rows*cols);
    
    for(i=0; i < rows; ++i)              // rows
       for(j=0; j < cols; ++j)          // columns
          out[px2(i,j, rows, cols)] =0; 
    
    for(i=0; i < rows; ++i)              // rows
    {
        for(j=0; j < cols; ++j)          // columns
        {
            sum = 0;                     // init to 0 before sum

            for(m=0; m < kRows; ++m)     // kernel rows
            {
                mm = kRows - 1 - m;      // row index of flipped kernel

                for(n=0; n < kCols; ++n) // kernel columns
                {
                    nn = kCols - 1 - n;  // column index of flipped kernel

                    // index of input signal, used for checking boundary
                    ii = i + (m - kCenterY);
                    jj = j + (n - kCenterX);

                    // ignore input samples which are out of bound
                    if( ii >= 0 && ii < rows && jj >= 0 && jj < cols)
                      //sum = sum+ in[ii+jj*rows]*kernel[mm][nn];
                        sum = sum+ in[px2(ii, jj, rows, cols)]*kernel[mm][nn];
                }
            }
            //out[i+j*rows] =sum;
            out[px2(i, j, rows, cols)]= sum;
        }
    }
    return out;
}


void conv(double* in, int rows, int cols, double **kernel,int kRows, int kCols, double** &out)
{
    int i,j,m,n;
    int ii,jj,mm,nn;
    double sum =0;
    // find center position of kernel (half of kernel size)
    int kCenterX = round((kCols+1)/2);
    int kCenterY = round((kRows+1)/2);
    out = double_matrix(rows, cols);
    
    for(i=0; i < rows; ++i)              // rows
       for(j=0; j < cols; ++j)          // columns
          out[i][j] =0; 
    
    for(i=0; i < rows; ++i)              // rows
    {
        for(j=0; j < cols; ++j)          // columns
        {
            sum = 0;                     // init to 0 before sum

            for(m=0; m < kRows; ++m)     // kernel rows
            {
                mm = kRows - 1 - m;      // row index of flipped kernel

                for(n=0; n < kCols; ++n) // kernel columns
                {
                    nn = kCols - 1 - n;  // column index of flipped kernel

                    // index of input signal, used for checking boundary
                    ii = i + (m - kCenterY);
                    jj = j + (n - kCenterX);

                    // ignore input samples which are out of bound
                    if( ii >= 0 && ii < rows && jj >= 0 && jj < cols && mm>=0 && mm<kRows && nn>=0 && nn<kCols)
                    out[i][j] += in[px2(ii,jj, rows, cols)] * kernel[mm][nn];
                }
            }
        }
    }
}
