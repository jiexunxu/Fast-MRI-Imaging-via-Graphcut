#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "math.h"

/* Matlab calling syntax:
 * [unaryTerms, binaryTerms, binaryTermsCounter]=fusionMoveUnaryAndBinaryTermsMex(...
 * unaryTerms(:), size(unaryTerms, 1), binaryTerms(:), binaryTermsCounter, size(binaryTerms, 1), 
 * nbrList(:)-1, size(nbrList, 1), real(x), imag(x), real(xCand), imag(xCand), labelMap, ...
 * power, threshold);
 */
double EsWeight;
double* binaryTerms;

double* unaryTerms;

double* labelMap;


inline double singlePairEs(double power, double Es_threshold, double x1real, double x1imag, double x2real, double x2imag)
{
    double difreal=x1real-x2real;double difimag=x1imag-x2imag;
    double E=pow(sqrt(difreal*difreal+difimag*difimag), power);
//  double E=pow(pow(abs(difreal), power)+pow(abs(difimag), power), 1.0/power);
    if(E<Es_threshold)
        return E*EsWeight;
    return Es_threshold*EsWeight;
}

void extractMaxSupport(const mxArray* supportPtrs, double* vxlSptCounts, int v1, int &spt_max, double &wmax_real, double &wmax_imag)
{
    int spt1=(int)vxlSptCounts[v1];
    if(spt1<=0)
        return;
    
    const mxArray* cellPtr=mxGetCell(supportPtrs, v1);
    double* v1_spt_ptr=mxGetPr(cellPtr);
    double wreal, wimag, wabs, wabs_max=0;
    for(int i=0;i<spt1;i++)
    {
        wreal=v1_spt_ptr[i*3+1];
        wimag=v1_spt_ptr[i*3+2];
        wabs=wreal*wreal+wimag*wimag;
        if(wabs>wabs_max)
        {
            wabs_max=wabs;
            spt_max=(int)v1_spt_ptr[i*3];
            wmax_real=wreal;
            wmax_imag=wimag;
        }
    }
}

void processNbrPairs(double power, double Es_threshold, double* x0real, double* x0imag, 
    mwSize unaryTermsCount, mwSize binaryTermsCount, const mxArray* supportPtrs, double* vxlSptCounts, 
    int v1, int v2, double* &unaryTermsOut, double* &binaryTermsOut, mwSize &binaryTermsIdx, double termWeight)
{    
    double x1real=x0real[v1];double x1imag=x0imag[v1];
    double x2real=x0real[v2];double x2imag=x0imag[v2];
    int s1=(int)vxlSpts[v1];int s2=(int)vxlSpts[v2];
    double w1real=vxlWeightsReal[v1];double w1imag=vxlWeightsImag[v1];
    double w2real=vxlWeightsReal[v2];double w2imag=vxlWeightsImag[v2];
    if(debugPrint)
    {
        printf("===Processing vertices %d and %d===\n", s1, s2);
        printf("s1=%d, w1real=%g, w1imag=%g\n", s1, w1real, w1imag);
        printf("s2=%d, w2real=%g, w2imag=%g\n", s2, w2real, w2imag);
    }
 //   int s1=-1;int s2=-1;
 //   double w1real, w1imag, w2real, w2imag;
 //   extractMaxSupport(supportPtrs, vxlSptCounts, v1, s1, w1real, w1imag);
 //   extractMaxSupport(supportPtrs, vxlSptCounts, v2, s2, w2real, w2imag);
    if((s1==-1)&&(s2!=-1))
    {
        unaryTermsOut[s2-1]+=termWeight*singlePairEs(power, Es_threshold, x1real, x1imag, x2real, x2imag);
        unaryTermsOut[s2-1+unaryTermsCount]+=termWeight*singlePairEs(power, Es_threshold, x1real, x1imag, x2real+w2real, x2imag+w2imag);
        if(debugPrint)
            printf("s1==[], %d %g %g\n", s2, unaryTermsOut[s2-1], unaryTermsOut[s2-1+unaryTermsCount]);
    }else if((s1!=-1)&&(s2==-1))
    {
        unaryTermsOut[s1-1]+=termWeight*singlePairEs(power, Es_threshold, x1real, x1imag, x2real, x2imag);
        unaryTermsOut[s1-1+unaryTermsCount]+=termWeight*singlePairEs(power, Es_threshold, x1real+w1real, x1imag+w1imag, x2real, x2imag);
        if(debugPrint)
            printf("s2==[], %d %g %g\n", s1, unaryTermsOut[s1-1], unaryTermsOut[s1-1+unaryTermsCount]);
    }else if((s1!=-1)&&(s2!=-1))
    {
        if(s1==s2)
        {
            unaryTermsOut[s1-1]+=termWeight*singlePairEs(power, Es_threshold, x1real, x1imag, x2real, x2imag);
            unaryTermsOut[s1-1+unaryTermsCount]+=termWeight*singlePairEs(power, Es_threshold, x1real+w1real, x1imag+w1imag, x2real+w2real, x2imag+w2imag);
            if(debugPrint)
                printf("s1==s2!=[], %d %g %g\n", s1, unaryTermsOut[s1-1], unaryTermsOut[s1-1+unaryTermsCount]);

        }else{
            binaryTermsOut[binaryTermsIdx-1]=s1;
            binaryTermsOut[binaryTermsIdx-1+binaryTermsCount]=s2;
            binaryTermsOut[binaryTermsIdx-1+binaryTermsCount*2]=termWeight*singlePairEs(power, Es_threshold, x1real, x1imag, x2real, x2imag);
            binaryTermsOut[binaryTermsIdx-1+binaryTermsCount*3]=termWeight*singlePairEs(power, Es_threshold, x1real+w1real, x1imag+w1imag, x2real, x2imag);
            binaryTermsOut[binaryTermsIdx-1+binaryTermsCount*4]=termWeight*singlePairEs(power, Es_threshold, x1real, x1imag, x2real+w2real, x2imag+w2imag);
            binaryTermsOut[binaryTermsIdx-1+binaryTermsCount*5]=termWeight*singlePairEs(power, Es_threshold, x1real+w1real, x1imag+w1imag, x2real+w2real, x2imag+w2imag);   
            if(debugPrint)
                printf("s1!=s2!=[], %d %d %g %g %g %g\n", s1, s2, binaryTermsOut[binaryTermsIdx-1+binaryTermsCount*2], 
                binaryTermsOut[binaryTermsIdx-1+binaryTermsCount*3], binaryTermsOut[binaryTermsIdx-1+binaryTermsCount*4], binaryTermsOut[binaryTermsIdx-1+binaryTermsCount*5]);

            binaryTermsIdx++;
         }
    }
}

void computeNbrTerms(double power, double Es_threshold, double* x0real, double* x0imag, 
    mwSize nbrListCount, double* nbrList, mwSize unaryTermsCount, mwSize binaryTermsCount, const mxArray* supportPtrs, 
    double* vxlSptCounts, double* nbrTermWeights, double* &unaryTermsOut, double* &binaryTermsOut, mwSize &binaryTermsIdx)
{
    int v1, v2;
    for(int i=0;i<nbrListCount;i++)
    {
        v1=(int)nbrList[i];v2=(int)nbrList[i+nbrListCount];
        double termWeight=nbrTermWeights[i];
        processNbrPairs(power, Es_threshold, x0real, x0imag, unaryTermsCount, binaryTermsCount, 
            supportPtrs, vxlSptCounts, v1, v2, unaryTermsOut, binaryTermsOut, binaryTermsIdx, termWeight);        
    }
}

void printArray(double* array, int dim)
{
    for(int i=0;i<dim;i++)
        printf("%g ", array[i]);
    printf("\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize unaryTermsCount;
    double* unaryTerms;
    mwSize binaryTermsIdx;
    mwSize binaryTermsCount;
    double* binaryTerms;
    mwSize nbrListCount; 
    double* nbrList; // nbrList[i] and nbrList[i+nbrListCount] are neighboring voxels
    double* x0real;
    double* x0imag;
    mwSize numOfDims;
    double* vxlSptCounts;
    const mxArray* supportPtrs;
    double Es_threshold;
    double power;
    double* nbrTermWeights;
    
    double* unaryTermsOut;
    double* binaryTermsOut;
    
    unaryTermsCount=(int)mxGetScalar(prhs[0]);    
    unaryTerms=mxGetPr(prhs[1]);
    binaryTermsIdx=(int)mxGetScalar(prhs[2]);
    binaryTermsCount=(int)mxGetScalar(prhs[3]);
    binaryTerms=mxGetPr(prhs[4]);
    nbrListCount=(int)mxGetScalar(prhs[5]);
    nbrList=mxGetPr(prhs[6]);
    x0real=mxGetPr(prhs[7]);
    x0imag=mxGetPr(prhs[8]);
    numOfDims=(int)mxGetScalar(prhs[9]);
    vxlSpts=mxGetPr(prhs[10]);
    vxlWeightsReal=mxGetPr(prhs[11]);
    vxlWeightsImag=mxGetPr(prhs[12]);
   // vxlSptCounts=mxGetPr(prhs[10]);
   // supportPtrs=prhs[11];
    Es_threshold=mxGetScalar(prhs[13]);
    power=mxGetScalar(prhs[14]);
    EsWeight=mxGetScalar(prhs[15]);
    nbrTermWeights=mxGetPr(prhs[16]);
    
    plhs[0]=mxCreateDoubleMatrix(unaryTermsCount*2, 1, mxREAL);
    unaryTermsOut=mxGetPr(plhs[0]);
    for(int i=0;i<unaryTermsCount*2;i++)
        unaryTermsOut[i]=unaryTerms[i];
    
    plhs[1]=mxCreateDoubleMatrix(binaryTermsCount*6, 1, mxREAL);
    binaryTermsOut=mxGetPr(plhs[1]);
    for(int i=0;i<binaryTermsIdx;i++)
        for(int j=0;j<6;j++)
            binaryTermsOut[i+binaryTermsCount*j]=binaryTerms[i+binaryTermsCount*j];

    computeNbrTerms(power, Es_threshold, x0real, x0imag, nbrListCount, nbrList, 
        unaryTermsCount, binaryTermsCount, supportPtrs, vxlSptCounts, nbrTermWeights,
        unaryTermsOut, binaryTermsOut, binaryTermsIdx);
    
    plhs[2]=mxCreateDoubleScalar(binaryTermsIdx);    
}