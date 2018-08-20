#include "mex.h"
#include "math.h"

double *data;
double *filteredData;
double *gaussianGrid;
mwSize numRows, numCols, numSlices;
mwSize windowSize, w;
double sigmaD, sigmaR;
    
void init()
{
    int i, j, k;    
    double distsqr, exponent;
    
    w=windowSize*2+1;
    gaussianGrid=mxGetPr(mxCreateDoubleMatrix(w*w*w, 1, mxREAL));   
    for(i=0;i<w;i++){
        for(j=0;j<w;j++){
            for(k=0;k<w;k++){
                distsqr=(i-windowSize)*(i-windowSize)+(j-windowSize)*(j-windowSize)+(k-windowSize)*(k-windowSize);
                exponent=-distsqr/(2*sigmaD*sigmaD);
                gaussianGrid[i+j*w+k*w*w]=exp(exponent);
            }
        }
    }
}

void bilateralFiltering3DSingleVoxel(int x, int y, int z)
{
    int i, j, k, idx;   
    double dataCtr, intExponent, localRegion, intensityWeights, intWtxLocalSum, intWtSum;
    dataCtr=data[z*numRows*numCols+y*numRows+x]; 
    intWtxLocalSum=0;
    intWtSum=0;
    for(i=0;i<w;i++){
        for(j=0;j<w;j++){
            for(k=0;k<w;k++){
                localRegion=data[x+i-windowSize+(y+j-windowSize)*numRows+(z+k-windowSize)*numRows*numCols];  
                intExponent=localRegion-dataCtr;
                intExponent=-intExponent*intExponent/(2*sigmaR*sigmaR);
                intensityWeights=exp(intExponent);
                intensityWeights=intensityWeights*gaussianGrid[i+j*w+k*w*w];
                intWtxLocalSum=intWtxLocalSum+localRegion*intensityWeights;
                intWtSum=intWtSum+intensityWeights;        
            }
        }
    }
    filteredData[z*numRows*numCols+y*numRows+x]=intWtxLocalSum/intWtSum;
}

void bilateralFiltering3D()
{
    int i, j, k;
    for(i=windowSize;i<numRows-windowSize;i++)
        for(j=windowSize;j<numCols-windowSize;j++)
            for(k=windowSize;k<numSlices-windowSize;k++)
                bilateralFiltering3DSingleVoxel(i, j, k);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    data=mxGetPr(prhs[0]);
    numRows=(mwSize)mxGetScalar(prhs[1]);
    numCols=(mwSize)mxGetScalar(prhs[2]);
    numSlices=(mwSize)mxGetScalar(prhs[3]);
    windowSize=(mwSize)mxGetScalar(prhs[4]);
    sigmaD=mxGetScalar(prhs[5]);
    sigmaR=mxGetScalar(prhs[6]);
    
    plhs[0]=mxCreateDoubleMatrix(numRows*numCols*numSlices, 1, mxREAL);
    filteredData=mxGetPr(plhs[0]);
    
    init();
    bilateralFiltering3D();
}