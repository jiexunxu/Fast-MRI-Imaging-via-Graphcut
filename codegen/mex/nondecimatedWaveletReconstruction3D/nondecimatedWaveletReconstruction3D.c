/*
 * nondecimatedWaveletReconstruction3D.c
 *
 * Code generation for function 'nondecimatedWaveletReconstruction3D'
 *
 * C source code generated on: Thu Sep 06 12:17:49 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nondecimatedWaveletReconstruction3D.h"

/* Type Definitions */
#ifndef struct_emxArray__common
#define struct_emxArray__common
typedef struct emxArray__common
{
    void *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray__common;
#endif
#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T
typedef struct emxArray_int32_T
{
    int32_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_int32_T;
#endif
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
typedef struct emxArray_real_T
{
    real_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_real_T;
#endif

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtDCInfo emlrtDCI = { 12, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo b_emlrtDCI = { 13, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo c_emlrtDCI = { 13, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo d_emlrtDCI = { 14, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo e_emlrtDCI = { 14, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo f_emlrtDCI = { 15, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo g_emlrtDCI = { 15, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo h_emlrtDCI = { 16, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo i_emlrtDCI = { 16, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo j_emlrtDCI = { 17, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo k_emlrtDCI = { 17, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo l_emlrtDCI = { 18, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo m_emlrtDCI = { 18, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo n_emlrtDCI = { 19, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo o_emlrtDCI = { 19, 34, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo p_emlrtDCI = { 20, 19, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo q_emlrtDCI = { 11, 27, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };
static emlrtDCInfo r_emlrtDCI = { 11, 27, "nondecimatedWaveletReconstruction3D", "X:/Projects/ApproximateNullVec/deployment/core/nondecimatedWaveletReconstruction3D.m", 1 };

/* Function Declarations */
static void b_conv(const creal_T A[183], const real_T B[8], creal_T C[190]);
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_creal_T *y);
static void b_emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void b_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B);
static void c_conv(const creal_T A[263], const real_T B[8], creal_T C[270]);
static void c_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_creal_T *ret);
static void c_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B);
static void conv(const emxArray_creal_T *A, emxArray_creal_T *C);
static void d_conv(const emxArray_creal_T *A, emxArray_creal_T *C);
static void d_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B);
static void e_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B);
static void emlrt_marshallIn(const mxArray *A, const char_T *identifier, emxArray_creal_T *y);
static const mxArray *emlrt_marshallOut(emxArray_creal_T *u);
static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize);
static void emxFree_creal_T(emxArray_creal_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void f_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B);
static void g_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B);
static void h_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B);
static void singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B);
static void squeeze(const emxArray_creal_T *a, emxArray_creal_T *b);

/* Function Definitions */

/*
 * 
 */
static void b_conv(const creal_T A[183], const real_T B[8], creal_T C[190])
{
    int32_T jC;
    int32_T jA1;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    for (jC = 0; jC < 190; jC++) {
        if (8 < jC + 2) {
            jA1 = jC;
        } else {
            jA1 = 7;
        }
        if (183 < jC + 1) {
            jA2 = 183;
        } else {
            jA2 = jC + 1;
        }
        s_re = 0.0;
        s_im = 0.0;
        for (jA1 -= 7; jA1 + 1 <= jA2; jA1++) {
            s_re += A[jA1].re * B[jC - jA1];
            s_im += A[jA1].im * B[jC - jA1];
        }
        C[jC].re = s_re;
        C[jC].im = s_im;
    }
}

static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_creal_T *y)
{
    c_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void b_emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
    emxArray_creal_T *emxArray;
    int32_T loop_ub;
    int32_T i;
    *pEmxArray = (emxArray_creal_T *)malloc(sizeof(emxArray_creal_T));
    if (doPush) {
        emlrtPushHeapReferenceStack((void *)pEmxArray, (void (*)(void *, boolean_T))emxFree_creal_T);
    }
    emxArray = *pEmxArray;
    emxArray->data = (creal_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = TRUE;
    loop_ub = numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        emxArray->size[i] = 0;
    }
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void b_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i5;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r6;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv1[3];
    int32_T b_loop_ub;
    int32_T i6;
    int32_T c_loop_ub;
    int32_T i7;
    emxArray_creal_T c_v;
    uint32_T d_loop_ub;
    uint32_T k;
    creal_T b[183];
    creal_T d_v[190];
    static const real_T dv2[8] = { 0.16290171402561984, 0.5054728575456503, 0.44610006912318972, -0.019787513117909966, -0.13225358368436993, 0.021808150237390023, 0.023251800535559971, -0.0074934946651300143 };
    creal_T b_b[263];
    creal_T e_v[270];
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'nondecimatedWaveletReconstruction3D:25' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:26' q=length(w1); */
    /* 'nondecimatedWaveletReconstruction3D:27' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    i5 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 270;
    B->size[1] = 190;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i5, (int32_T)sizeof(creal_T));
    loop_ub = 51300 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i5 = 0; i5 <= loop_ub; i5++) {
        B->data[i5].re = 0.0;
        B->data[i5].im = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:28' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r6, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:29' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletReconstruction3D:30' v=conv(squeeze(A(i, j, :)), w1); */
            i5 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i5, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i5 = 0; i5 <= loop_ub; i5++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i5] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i5];
            }
            squeeze(b_A, v);
            i5 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i5, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i5 = 0; i5 <= loop_ub; i5++) {
                b_v->data[i5] = v->data[i5];
            }
            d_conv(b_v, v);
            /* 'nondecimatedWaveletReconstruction3D:31' B(i, j, :)=v; */
            i5 = r6->size[0];
            r6->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r6, i5, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (i5 = 0; i5 <= loop_ub; i5++) {
                r6->data[i5] = i5;
            }
            iv1[0] = 1;
            iv1[1] = 1;
            iv1[2] = r6->size[0];
            loop_ub = iv1[2] - 1;
            for (i5 = 0; i5 <= loop_ub; i5++) {
                b_loop_ub = iv1[1] - 1;
                for (i6 = 0; i6 <= b_loop_ub; i6++) {
                    c_loop_ub = iv1[0] - 1;
                    for (i7 = 0; i7 <= c_loop_ub; i7++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv1;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r6->data[i5]] = c_v.data[(i7 + c_v.size[0] * i6) + c_v.size[0] * c_v.size[1] * i5];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r6);
    emxFree_creal_T(&v);
    /* 'nondecimatedWaveletReconstruction3D:35' for i=1:m */
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:36' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:37' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (loop_ub = 0; loop_ub < 183; loop_ub++) {
                b[loop_ub] = B->data[(i + B->size[0] * loop_ub) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv2, d_v);
            /* 'nondecimatedWaveletReconstruction3D:38' B(i, :, k)=v; */
            loop_ub = (int32_T)k;
            for (i5 = 0; i5 < 190; i5++) {
                B->data[(i + B->size[0] * i5) + B->size[0] * B->size[1] * (loop_ub - 1)] = d_v[i5];
            }
        }
    }
    /* 'nondecimatedWaveletReconstruction3D:42' for j=1:n+q-1 */
    for (j = 0; j < 190; j++) {
        /* 'nondecimatedWaveletReconstruction3D:43' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:44' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (loop_ub = 0; loop_ub < 263; loop_ub++) {
                b_b[loop_ub] = B->data[(loop_ub + B->size[0] * j) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            c_conv(b_b, dv2, e_v);
            /* 'nondecimatedWaveletReconstruction3D:45' B(:, j, k)=v; */
            loop_ub = (int32_T)k;
            for (i5 = 0; i5 < 270; i5++) {
                B->data[(i5 + B->size[0] * j) + B->size[0] * B->size[1] * (loop_ub - 1)] = e_v[i5];
            }
        }
    }
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * 
 */
static void c_conv(const creal_T A[263], const real_T B[8], creal_T C[270])
{
    int32_T jC;
    int32_T jA1;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    for (jC = 0; jC < 270; jC++) {
        if (8 < jC + 2) {
            jA1 = jC;
        } else {
            jA1 = 7;
        }
        if (263 < jC + 1) {
            jA2 = 263;
        } else {
            jA2 = jC + 1;
        }
        s_re = 0.0;
        s_im = 0.0;
        for (jA1 -= 7; jA1 + 1 <= jA2; jA1++) {
            s_re += A[jA1].re * B[jC - jA1];
            s_im += A[jA1].im * B[jC - jA1];
        }
        C[jC].re = s_re;
        C[jC].im = s_im;
    }
}

static void c_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_creal_T *ret)
{
    int32_T i;
    static const int16_T iv8[3] = { 263, 183, 2104 };
    int32_T iv9[3];
    static const boolean_T bv0[3] = { FALSE, FALSE, TRUE };
    boolean_T bv1[3];
    for (i = 0; i < 3; i++) {
        iv9[i] = iv8[i];
        bv1[i] = bv0[i];
    }
    emlrtCheckVsBuiltInR2011a(msgId, src, "double", TRUE, 3U, iv9, bv1, ret->size);
    i = ret->size[0] * ret->size[1] * ret->size[2];
    ret->size[0] = 263;
    ret->size[1] = 183;
    emxEnsureCapacity((emxArray__common *)ret, i, (int32_T)sizeof(creal_T));
    emlrtImportArrayR2008b(src, ret->data, 8);
    emlrtDestroyArray(&src);
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void c_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i8;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r7;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv2[3];
    int32_T b_loop_ub;
    int32_T i9;
    int32_T c_loop_ub;
    int32_T i10;
    emxArray_creal_T c_v;
    uint32_T d_loop_ub;
    uint32_T k;
    creal_T b[183];
    creal_T d_v[190];
    static const real_T dv4[8] = { -0.0074934946651300143, -0.023251800535559971, 0.021808150237390023, 0.13225358368436993, -0.019787513117909966, -0.44610006912318972, 0.5054728575456503, -0.16290171402561984 };
    creal_T b_b[263];
    creal_T e_v[270];
    static const real_T dv5[8] = { 0.16290171402561984, 0.5054728575456503, 0.44610006912318972, -0.019787513117909966, -0.13225358368436993, 0.021808150237390023, 0.023251800535559971, -0.0074934946651300143 };
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'nondecimatedWaveletReconstruction3D:25' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:26' q=length(w1); */
    /* 'nondecimatedWaveletReconstruction3D:27' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    i8 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 270;
    B->size[1] = 190;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i8, (int32_T)sizeof(creal_T));
    loop_ub = 51300 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        B->data[i8].re = 0.0;
        B->data[i8].im = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:28' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r7, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:29' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletReconstruction3D:30' v=conv(squeeze(A(i, j, :)), w1); */
            i8 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i8, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i8] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i8];
            }
            squeeze(b_A, v);
            i8 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i8, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                b_v->data[i8] = v->data[i8];
            }
            conv(b_v, v);
            /* 'nondecimatedWaveletReconstruction3D:31' B(i, j, :)=v; */
            i8 = r7->size[0];
            r7->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r7, i8, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                r7->data[i8] = i8;
            }
            iv2[0] = 1;
            iv2[1] = 1;
            iv2[2] = r7->size[0];
            loop_ub = iv2[2] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                b_loop_ub = iv2[1] - 1;
                for (i9 = 0; i9 <= b_loop_ub; i9++) {
                    c_loop_ub = iv2[0] - 1;
                    for (i10 = 0; i10 <= c_loop_ub; i10++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv2;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r7->data[i8]] = c_v.data[(i10 + c_v.size[0] * i9) + c_v.size[0] * c_v.size[1] * i8];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r7);
    emxFree_creal_T(&v);
    /* 'nondecimatedWaveletReconstruction3D:35' for i=1:m */
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:36' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:37' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (loop_ub = 0; loop_ub < 183; loop_ub++) {
                b[loop_ub] = B->data[(i + B->size[0] * loop_ub) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv4, d_v);
            /* 'nondecimatedWaveletReconstruction3D:38' B(i, :, k)=v; */
            loop_ub = (int32_T)k;
            for (i8 = 0; i8 < 190; i8++) {
                B->data[(i + B->size[0] * i8) + B->size[0] * B->size[1] * (loop_ub - 1)] = d_v[i8];
            }
        }
    }
    /* 'nondecimatedWaveletReconstruction3D:42' for j=1:n+q-1 */
    for (j = 0; j < 190; j++) {
        /* 'nondecimatedWaveletReconstruction3D:43' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:44' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (loop_ub = 0; loop_ub < 263; loop_ub++) {
                b_b[loop_ub] = B->data[(loop_ub + B->size[0] * j) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            c_conv(b_b, dv5, e_v);
            /* 'nondecimatedWaveletReconstruction3D:45' B(:, j, k)=v; */
            loop_ub = (int32_T)k;
            for (i8 = 0; i8 < 270; i8++) {
                B->data[(i8 + B->size[0] * j) + B->size[0] * B->size[1] * (loop_ub - 1)] = e_v[i8];
            }
        }
    }
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * 
 */
static void conv(const emxArray_creal_T *A, emxArray_creal_T *C)
{
    int32_T nA;
    int32_T nApnB;
    int32_T jA1;
    int32_T jC;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    static const real_T dv1[8] = { 0.16290171402561984, 0.5054728575456503, 0.44610006912318972, -0.019787513117909966, -0.13225358368436993, 0.021808150237390023, 0.023251800535559971, -0.0074934946651300143 };
    nA = A->size[0];
    nApnB = nA + 7;
    if (nA == 0) {
        nApnB++;
    }
    jA1 = C->size[0];
    C->size[0] = nApnB;
    emxEnsureCapacity((emxArray__common *)C, jA1, (int32_T)sizeof(creal_T));
    if ((A->size[0] == 0) || (nApnB == 0)) {
        jC = C->size[0];
        jA1 = C->size[0];
        C->size[0] = jC;
        emxEnsureCapacity((emxArray__common *)C, jA1, (int32_T)sizeof(creal_T));
        jC--;
        for (jA1 = 0; jA1 <= jC; jA1++) {
            C->data[jA1].re = 0.0;
            C->data[jA1].im = 0.0;
        }
    } else {
        for (jC = 1; jC <= nApnB; jC++) {
            if (8 < jC + 1) {
                jA1 = jC - 7;
            } else {
                jA1 = 1;
            }
            if (nA < jC) {
                jA2 = nA;
            } else {
                jA2 = jC;
            }
            s_re = 0.0;
            s_im = 0.0;
            while (jA1 <= jA2) {
                s_re += A->data[jA1 - 1].re * dv1[jC - jA1];
                s_im += A->data[jA1 - 1].im * dv1[jC - jA1];
                jA1++;
            }
            C->data[jC - 1].re = s_re;
            C->data[jC - 1].im = s_im;
        }
    }
}

/*
 * 
 */
static void d_conv(const emxArray_creal_T *A, emxArray_creal_T *C)
{
    int32_T nA;
    int32_T nApnB;
    int32_T jA1;
    int32_T jC;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    static const real_T dv3[8] = { -0.0074934946651300143, -0.023251800535559971, 0.021808150237390023, 0.13225358368436993, -0.019787513117909966, -0.44610006912318972, 0.5054728575456503, -0.16290171402561984 };
    nA = A->size[0];
    nApnB = nA + 7;
    if (nA == 0) {
        nApnB++;
    }
    jA1 = C->size[0];
    C->size[0] = nApnB;
    emxEnsureCapacity((emxArray__common *)C, jA1, (int32_T)sizeof(creal_T));
    if ((A->size[0] == 0) || (nApnB == 0)) {
        jC = C->size[0];
        jA1 = C->size[0];
        C->size[0] = jC;
        emxEnsureCapacity((emxArray__common *)C, jA1, (int32_T)sizeof(creal_T));
        jC--;
        for (jA1 = 0; jA1 <= jC; jA1++) {
            C->data[jA1].re = 0.0;
            C->data[jA1].im = 0.0;
        }
    } else {
        for (jC = 1; jC <= nApnB; jC++) {
            if (8 < jC + 1) {
                jA1 = jC - 7;
            } else {
                jA1 = 1;
            }
            if (nA < jC) {
                jA2 = nA;
            } else {
                jA2 = jC;
            }
            s_re = 0.0;
            s_im = 0.0;
            while (jA1 <= jA2) {
                s_re += A->data[jA1 - 1].re * dv3[jC - jA1];
                s_im += A->data[jA1 - 1].im * dv3[jC - jA1];
                jA1++;
            }
            C->data[jC - 1].re = s_re;
            C->data[jC - 1].im = s_im;
        }
    }
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void d_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i11;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r8;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv3[3];
    int32_T b_loop_ub;
    int32_T i12;
    int32_T c_loop_ub;
    int32_T i13;
    emxArray_creal_T c_v;
    uint32_T d_loop_ub;
    uint32_T k;
    creal_T b[183];
    creal_T d_v[190];
    static const real_T dv6[8] = { -0.0074934946651300143, -0.023251800535559971, 0.021808150237390023, 0.13225358368436993, -0.019787513117909966, -0.44610006912318972, 0.5054728575456503, -0.16290171402561984 };
    creal_T b_b[263];
    creal_T e_v[270];
    static const real_T dv7[8] = { 0.16290171402561984, 0.5054728575456503, 0.44610006912318972, -0.019787513117909966, -0.13225358368436993, 0.021808150237390023, 0.023251800535559971, -0.0074934946651300143 };
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'nondecimatedWaveletReconstruction3D:25' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:26' q=length(w1); */
    /* 'nondecimatedWaveletReconstruction3D:27' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    i11 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 270;
    B->size[1] = 190;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i11, (int32_T)sizeof(creal_T));
    loop_ub = 51300 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i11 = 0; i11 <= loop_ub; i11++) {
        B->data[i11].re = 0.0;
        B->data[i11].im = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:28' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r8, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:29' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletReconstruction3D:30' v=conv(squeeze(A(i, j, :)), w1); */
            i11 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i11, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i11 = 0; i11 <= loop_ub; i11++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i11] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i11];
            }
            squeeze(b_A, v);
            i11 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i11, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i11 = 0; i11 <= loop_ub; i11++) {
                b_v->data[i11] = v->data[i11];
            }
            d_conv(b_v, v);
            /* 'nondecimatedWaveletReconstruction3D:31' B(i, j, :)=v; */
            i11 = r8->size[0];
            r8->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r8, i11, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (i11 = 0; i11 <= loop_ub; i11++) {
                r8->data[i11] = i11;
            }
            iv3[0] = 1;
            iv3[1] = 1;
            iv3[2] = r8->size[0];
            loop_ub = iv3[2] - 1;
            for (i11 = 0; i11 <= loop_ub; i11++) {
                b_loop_ub = iv3[1] - 1;
                for (i12 = 0; i12 <= b_loop_ub; i12++) {
                    c_loop_ub = iv3[0] - 1;
                    for (i13 = 0; i13 <= c_loop_ub; i13++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv3;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r8->data[i11]] = c_v.data[(i13 + c_v.size[0] * i12) + c_v.size[0] * c_v.size[1] * i11];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r8);
    emxFree_creal_T(&v);
    /* 'nondecimatedWaveletReconstruction3D:35' for i=1:m */
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:36' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:37' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (loop_ub = 0; loop_ub < 183; loop_ub++) {
                b[loop_ub] = B->data[(i + B->size[0] * loop_ub) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv6, d_v);
            /* 'nondecimatedWaveletReconstruction3D:38' B(i, :, k)=v; */
            loop_ub = (int32_T)k;
            for (i11 = 0; i11 < 190; i11++) {
                B->data[(i + B->size[0] * i11) + B->size[0] * B->size[1] * (loop_ub - 1)] = d_v[i11];
            }
        }
    }
    /* 'nondecimatedWaveletReconstruction3D:42' for j=1:n+q-1 */
    for (j = 0; j < 190; j++) {
        /* 'nondecimatedWaveletReconstruction3D:43' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:44' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (loop_ub = 0; loop_ub < 263; loop_ub++) {
                b_b[loop_ub] = B->data[(loop_ub + B->size[0] * j) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            c_conv(b_b, dv7, e_v);
            /* 'nondecimatedWaveletReconstruction3D:45' B(:, j, k)=v; */
            loop_ub = (int32_T)k;
            for (i11 = 0; i11 < 270; i11++) {
                B->data[(i11 + B->size[0] * j) + B->size[0] * B->size[1] * (loop_ub - 1)] = e_v[i11];
            }
        }
    }
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void e_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i14;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r9;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv4[3];
    int32_T b_loop_ub;
    int32_T i15;
    int32_T c_loop_ub;
    int32_T i16;
    emxArray_creal_T c_v;
    uint32_T d_loop_ub;
    uint32_T k;
    creal_T b[183];
    creal_T d_v[190];
    static const real_T dv8[8] = { 0.16290171402561984, 0.5054728575456503, 0.44610006912318972, -0.019787513117909966, -0.13225358368436993, 0.021808150237390023, 0.023251800535559971, -0.0074934946651300143 };
    creal_T b_b[263];
    creal_T e_v[270];
    static const real_T dv9[8] = { -0.0074934946651300143, -0.023251800535559971, 0.021808150237390023, 0.13225358368436993, -0.019787513117909966, -0.44610006912318972, 0.5054728575456503, -0.16290171402561984 };
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'nondecimatedWaveletReconstruction3D:25' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:26' q=length(w1); */
    /* 'nondecimatedWaveletReconstruction3D:27' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    i14 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 270;
    B->size[1] = 190;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i14, (int32_T)sizeof(creal_T));
    loop_ub = 51300 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i14 = 0; i14 <= loop_ub; i14++) {
        B->data[i14].re = 0.0;
        B->data[i14].im = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:28' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r9, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:29' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletReconstruction3D:30' v=conv(squeeze(A(i, j, :)), w1); */
            i14 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i14, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i14 = 0; i14 <= loop_ub; i14++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i14] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i14];
            }
            squeeze(b_A, v);
            i14 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i14, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i14 = 0; i14 <= loop_ub; i14++) {
                b_v->data[i14] = v->data[i14];
            }
            conv(b_v, v);
            /* 'nondecimatedWaveletReconstruction3D:31' B(i, j, :)=v; */
            i14 = r9->size[0];
            r9->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r9, i14, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (i14 = 0; i14 <= loop_ub; i14++) {
                r9->data[i14] = i14;
            }
            iv4[0] = 1;
            iv4[1] = 1;
            iv4[2] = r9->size[0];
            loop_ub = iv4[2] - 1;
            for (i14 = 0; i14 <= loop_ub; i14++) {
                b_loop_ub = iv4[1] - 1;
                for (i15 = 0; i15 <= b_loop_ub; i15++) {
                    c_loop_ub = iv4[0] - 1;
                    for (i16 = 0; i16 <= c_loop_ub; i16++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv4;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r9->data[i14]] = c_v.data[(i16 + c_v.size[0] * i15) + c_v.size[0] * c_v.size[1] * i14];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r9);
    emxFree_creal_T(&v);
    /* 'nondecimatedWaveletReconstruction3D:35' for i=1:m */
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:36' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:37' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (loop_ub = 0; loop_ub < 183; loop_ub++) {
                b[loop_ub] = B->data[(i + B->size[0] * loop_ub) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv8, d_v);
            /* 'nondecimatedWaveletReconstruction3D:38' B(i, :, k)=v; */
            loop_ub = (int32_T)k;
            for (i14 = 0; i14 < 190; i14++) {
                B->data[(i + B->size[0] * i14) + B->size[0] * B->size[1] * (loop_ub - 1)] = d_v[i14];
            }
        }
    }
    /* 'nondecimatedWaveletReconstruction3D:42' for j=1:n+q-1 */
    for (j = 0; j < 190; j++) {
        /* 'nondecimatedWaveletReconstruction3D:43' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:44' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (loop_ub = 0; loop_ub < 263; loop_ub++) {
                b_b[loop_ub] = B->data[(loop_ub + B->size[0] * j) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            c_conv(b_b, dv9, e_v);
            /* 'nondecimatedWaveletReconstruction3D:45' B(:, j, k)=v; */
            loop_ub = (int32_T)k;
            for (i14 = 0; i14 < 270; i14++) {
                B->data[(i14 + B->size[0] * j) + B->size[0] * B->size[1] * (loop_ub - 1)] = e_v[i14];
            }
        }
    }
    emlrtHeapReferenceStackLeaveFcn();
}

static void emlrt_marshallIn(const mxArray *A, const char_T *identifier, emxArray_creal_T *y)
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    b_emlrt_marshallIn(emlrtAlias(A), &thisId, y);
    emlrtDestroyArray(&A);
}

static const mxArray *emlrt_marshallOut(emxArray_creal_T *u)
{
    const mxArray *y;
    const mxArray *m0;
    y = NULL;
    m0 = mxCreateNumericArray(3, u->size, mxDOUBLE_CLASS, mxCOMPLEX);
    emlrtExportNumericArrayR2008b(m0, (void *)u->data, 8);
    emlrtAssign(&y, m0);
    return y;
}

static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize)
{
    int32_T newNumel;
    int32_T loop_ub;
    int32_T i;
    void *newData;
    newNumel = 1;
    loop_ub = emxArray->numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        newNumel *= emxArray->size[i];
    }
    if (newNumel > emxArray->allocatedSize) {
        loop_ub = emxArray->allocatedSize;
        if (loop_ub < 16) {
            loop_ub = 16;
        }
        while (loop_ub < newNumel) {
            loop_ub <<= 1;
        }
        newData = calloc((uint32_T)loop_ub, (uint32_T)elementSize);
        if (emxArray->data != NULL) {
            memcpy(newData, emxArray->data, (uint32_T)(elementSize * oldNumel));
            if (emxArray->canFreeData) {
                free(emxArray->data);
            }
        }
        emxArray->data = newData;
        emxArray->allocatedSize = loop_ub;
        emxArray->canFreeData = TRUE;
    }
}

static void emxFree_creal_T(emxArray_creal_T **pEmxArray)
{
    if (*pEmxArray != (emxArray_creal_T *)NULL) {
        if ((*pEmxArray)->canFreeData) {
            free((void *)(*pEmxArray)->data);
        }
        free((void *)(*pEmxArray)->size);
        free((void *)*pEmxArray);
        *pEmxArray = (emxArray_creal_T *)NULL;
    }
}

static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
    if (*pEmxArray != (emxArray_int32_T *)NULL) {
        if ((*pEmxArray)->canFreeData) {
            free((void *)(*pEmxArray)->data);
        }
        free((void *)(*pEmxArray)->size);
        free((void *)*pEmxArray);
        *pEmxArray = (emxArray_int32_T *)NULL;
    }
}

static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
    if (*pEmxArray != (emxArray_real_T *)NULL) {
        if ((*pEmxArray)->canFreeData) {
            free((void *)(*pEmxArray)->data);
        }
        free((void *)(*pEmxArray)->size);
        free((void *)*pEmxArray);
        *pEmxArray = (emxArray_real_T *)NULL;
    }
}

static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
    emxArray_creal_T *emxArray;
    int32_T loop_ub;
    int32_T i;
    *pEmxArray = (emxArray_creal_T *)malloc(sizeof(emxArray_creal_T));
    if (doPush) {
        emlrtPushHeapReferenceStack((void *)pEmxArray, (void (*)(void *, boolean_T))emxFree_creal_T);
    }
    emxArray = *pEmxArray;
    emxArray->data = (creal_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = TRUE;
    loop_ub = numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        emxArray->size[i] = 0;
    }
}

static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
    emxArray_int32_T *emxArray;
    int32_T loop_ub;
    int32_T i;
    *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
    if (doPush) {
        emlrtPushHeapReferenceStack((void *)pEmxArray, (void (*)(void *, boolean_T))emxFree_int32_T);
    }
    emxArray = *pEmxArray;
    emxArray->data = (int32_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = TRUE;
    loop_ub = numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        emxArray->size[i] = 0;
    }
}

static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
    emxArray_real_T *emxArray;
    int32_T loop_ub;
    int32_T i;
    *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
    if (doPush) {
        emlrtPushHeapReferenceStack((void *)pEmxArray, (void (*)(void *, boolean_T))emxFree_real_T);
    }
    emxArray = *pEmxArray;
    emxArray->data = (real_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = TRUE;
    loop_ub = numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        emxArray->size[i] = 0;
    }
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void f_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i17;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r10;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv5[3];
    int32_T b_loop_ub;
    int32_T i18;
    int32_T c_loop_ub;
    int32_T i19;
    emxArray_creal_T c_v;
    uint32_T d_loop_ub;
    uint32_T k;
    creal_T b[183];
    creal_T d_v[190];
    static const real_T dv10[8] = { 0.16290171402561984, 0.5054728575456503, 0.44610006912318972, -0.019787513117909966, -0.13225358368436993, 0.021808150237390023, 0.023251800535559971, -0.0074934946651300143 };
    creal_T b_b[263];
    creal_T e_v[270];
    static const real_T dv11[8] = { -0.0074934946651300143, -0.023251800535559971, 0.021808150237390023, 0.13225358368436993, -0.019787513117909966, -0.44610006912318972, 0.5054728575456503, -0.16290171402561984 };
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'nondecimatedWaveletReconstruction3D:25' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:26' q=length(w1); */
    /* 'nondecimatedWaveletReconstruction3D:27' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    i17 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 270;
    B->size[1] = 190;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i17, (int32_T)sizeof(creal_T));
    loop_ub = 51300 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i17 = 0; i17 <= loop_ub; i17++) {
        B->data[i17].re = 0.0;
        B->data[i17].im = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:28' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r10, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:29' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletReconstruction3D:30' v=conv(squeeze(A(i, j, :)), w1); */
            i17 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i17, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i17 = 0; i17 <= loop_ub; i17++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i17] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i17];
            }
            squeeze(b_A, v);
            i17 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i17, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i17 = 0; i17 <= loop_ub; i17++) {
                b_v->data[i17] = v->data[i17];
            }
            d_conv(b_v, v);
            /* 'nondecimatedWaveletReconstruction3D:31' B(i, j, :)=v; */
            i17 = r10->size[0];
            r10->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r10, i17, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (i17 = 0; i17 <= loop_ub; i17++) {
                r10->data[i17] = i17;
            }
            iv5[0] = 1;
            iv5[1] = 1;
            iv5[2] = r10->size[0];
            loop_ub = iv5[2] - 1;
            for (i17 = 0; i17 <= loop_ub; i17++) {
                b_loop_ub = iv5[1] - 1;
                for (i18 = 0; i18 <= b_loop_ub; i18++) {
                    c_loop_ub = iv5[0] - 1;
                    for (i19 = 0; i19 <= c_loop_ub; i19++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv5;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r10->data[i17]] = c_v.data[(i19 + c_v.size[0] * i18) + c_v.size[0] * c_v.size[1] * i17];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r10);
    emxFree_creal_T(&v);
    /* 'nondecimatedWaveletReconstruction3D:35' for i=1:m */
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:36' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:37' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (loop_ub = 0; loop_ub < 183; loop_ub++) {
                b[loop_ub] = B->data[(i + B->size[0] * loop_ub) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv10, d_v);
            /* 'nondecimatedWaveletReconstruction3D:38' B(i, :, k)=v; */
            loop_ub = (int32_T)k;
            for (i17 = 0; i17 < 190; i17++) {
                B->data[(i + B->size[0] * i17) + B->size[0] * B->size[1] * (loop_ub - 1)] = d_v[i17];
            }
        }
    }
    /* 'nondecimatedWaveletReconstruction3D:42' for j=1:n+q-1 */
    for (j = 0; j < 190; j++) {
        /* 'nondecimatedWaveletReconstruction3D:43' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:44' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (loop_ub = 0; loop_ub < 263; loop_ub++) {
                b_b[loop_ub] = B->data[(loop_ub + B->size[0] * j) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            c_conv(b_b, dv11, e_v);
            /* 'nondecimatedWaveletReconstruction3D:45' B(:, j, k)=v; */
            loop_ub = (int32_T)k;
            for (i17 = 0; i17 < 270; i17++) {
                B->data[(i17 + B->size[0] * j) + B->size[0] * B->size[1] * (loop_ub - 1)] = e_v[i17];
            }
        }
    }
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void g_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i20;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r11;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv6[3];
    int32_T b_loop_ub;
    int32_T i21;
    int32_T c_loop_ub;
    int32_T i22;
    emxArray_creal_T c_v;
    uint32_T d_loop_ub;
    uint32_T k;
    creal_T b[183];
    creal_T d_v[190];
    static const real_T dv12[8] = { -0.0074934946651300143, -0.023251800535559971, 0.021808150237390023, 0.13225358368436993, -0.019787513117909966, -0.44610006912318972, 0.5054728575456503, -0.16290171402561984 };
    creal_T b_b[263];
    creal_T e_v[270];
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'nondecimatedWaveletReconstruction3D:25' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:26' q=length(w1); */
    /* 'nondecimatedWaveletReconstruction3D:27' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    i20 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 270;
    B->size[1] = 190;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i20, (int32_T)sizeof(creal_T));
    loop_ub = 51300 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i20 = 0; i20 <= loop_ub; i20++) {
        B->data[i20].re = 0.0;
        B->data[i20].im = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:28' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r11, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:29' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletReconstruction3D:30' v=conv(squeeze(A(i, j, :)), w1); */
            i20 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i20, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i20 = 0; i20 <= loop_ub; i20++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i20] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i20];
            }
            squeeze(b_A, v);
            i20 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i20, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i20 = 0; i20 <= loop_ub; i20++) {
                b_v->data[i20] = v->data[i20];
            }
            conv(b_v, v);
            /* 'nondecimatedWaveletReconstruction3D:31' B(i, j, :)=v; */
            i20 = r11->size[0];
            r11->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r11, i20, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (i20 = 0; i20 <= loop_ub; i20++) {
                r11->data[i20] = i20;
            }
            iv6[0] = 1;
            iv6[1] = 1;
            iv6[2] = r11->size[0];
            loop_ub = iv6[2] - 1;
            for (i20 = 0; i20 <= loop_ub; i20++) {
                b_loop_ub = iv6[1] - 1;
                for (i21 = 0; i21 <= b_loop_ub; i21++) {
                    c_loop_ub = iv6[0] - 1;
                    for (i22 = 0; i22 <= c_loop_ub; i22++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv6;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r11->data[i20]] = c_v.data[(i22 + c_v.size[0] * i21) + c_v.size[0] * c_v.size[1] * i20];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r11);
    emxFree_creal_T(&v);
    /* 'nondecimatedWaveletReconstruction3D:35' for i=1:m */
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:36' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:37' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (loop_ub = 0; loop_ub < 183; loop_ub++) {
                b[loop_ub] = B->data[(i + B->size[0] * loop_ub) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv12, d_v);
            /* 'nondecimatedWaveletReconstruction3D:38' B(i, :, k)=v; */
            loop_ub = (int32_T)k;
            for (i20 = 0; i20 < 190; i20++) {
                B->data[(i + B->size[0] * i20) + B->size[0] * B->size[1] * (loop_ub - 1)] = d_v[i20];
            }
        }
    }
    /* 'nondecimatedWaveletReconstruction3D:42' for j=1:n+q-1 */
    for (j = 0; j < 190; j++) {
        /* 'nondecimatedWaveletReconstruction3D:43' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:44' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (loop_ub = 0; loop_ub < 263; loop_ub++) {
                b_b[loop_ub] = B->data[(loop_ub + B->size[0] * j) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            c_conv(b_b, dv12, e_v);
            /* 'nondecimatedWaveletReconstruction3D:45' B(:, j, k)=v; */
            loop_ub = (int32_T)k;
            for (i20 = 0; i20 < 270; i20++) {
                B->data[(i20 + B->size[0] * j) + B->size[0] * B->size[1] * (loop_ub - 1)] = e_v[i20];
            }
        }
    }
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void h_singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i23;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r12;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv7[3];
    int32_T b_loop_ub;
    int32_T i24;
    int32_T c_loop_ub;
    int32_T i25;
    emxArray_creal_T c_v;
    uint32_T d_loop_ub;
    uint32_T k;
    creal_T b[183];
    creal_T d_v[190];
    static const real_T dv13[8] = { -0.0074934946651300143, -0.023251800535559971, 0.021808150237390023, 0.13225358368436993, -0.019787513117909966, -0.44610006912318972, 0.5054728575456503, -0.16290171402561984 };
    creal_T b_b[263];
    creal_T e_v[270];
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'nondecimatedWaveletReconstruction3D:25' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:26' q=length(w1); */
    /* 'nondecimatedWaveletReconstruction3D:27' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    i23 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 270;
    B->size[1] = 190;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i23, (int32_T)sizeof(creal_T));
    loop_ub = 51300 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i23 = 0; i23 <= loop_ub; i23++) {
        B->data[i23].re = 0.0;
        B->data[i23].im = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:28' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r12, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:29' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletReconstruction3D:30' v=conv(squeeze(A(i, j, :)), w1); */
            i23 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i23, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i23 = 0; i23 <= loop_ub; i23++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i23] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i23];
            }
            squeeze(b_A, v);
            i23 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i23, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i23 = 0; i23 <= loop_ub; i23++) {
                b_v->data[i23] = v->data[i23];
            }
            d_conv(b_v, v);
            /* 'nondecimatedWaveletReconstruction3D:31' B(i, j, :)=v; */
            i23 = r12->size[0];
            r12->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r12, i23, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (i23 = 0; i23 <= loop_ub; i23++) {
                r12->data[i23] = i23;
            }
            iv7[0] = 1;
            iv7[1] = 1;
            iv7[2] = r12->size[0];
            loop_ub = iv7[2] - 1;
            for (i23 = 0; i23 <= loop_ub; i23++) {
                b_loop_ub = iv7[1] - 1;
                for (i24 = 0; i24 <= b_loop_ub; i24++) {
                    c_loop_ub = iv7[0] - 1;
                    for (i25 = 0; i25 <= c_loop_ub; i25++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv7;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r12->data[i23]] = c_v.data[(i25 + c_v.size[0] * i24) + c_v.size[0] * c_v.size[1] * i23];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r12);
    emxFree_creal_T(&v);
    /* 'nondecimatedWaveletReconstruction3D:35' for i=1:m */
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:36' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:37' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (loop_ub = 0; loop_ub < 183; loop_ub++) {
                b[loop_ub] = B->data[(i + B->size[0] * loop_ub) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv13, d_v);
            /* 'nondecimatedWaveletReconstruction3D:38' B(i, :, k)=v; */
            loop_ub = (int32_T)k;
            for (i23 = 0; i23 < 190; i23++) {
                B->data[(i + B->size[0] * i23) + B->size[0] * B->size[1] * (loop_ub - 1)] = d_v[i23];
            }
        }
    }
    /* 'nondecimatedWaveletReconstruction3D:42' for j=1:n+q-1 */
    for (j = 0; j < 190; j++) {
        /* 'nondecimatedWaveletReconstruction3D:43' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:44' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (loop_ub = 0; loop_ub < 263; loop_ub++) {
                b_b[loop_ub] = B->data[(loop_ub + B->size[0] * j) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            c_conv(b_b, dv13, e_v);
            /* 'nondecimatedWaveletReconstruction3D:45' B(:, j, k)=v; */
            loop_ub = (int32_T)k;
            for (i23 = 0; i23 < 270; i23++) {
                B->data[(i23 + B->size[0] * j) + B->size[0] * B->size[1] * (loop_ub - 1)] = e_v[i23];
            }
        }
    }
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void singleBlockRecon(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i2;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r5;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv0[3];
    int32_T b_loop_ub;
    int32_T i3;
    int32_T c_loop_ub;
    int32_T i4;
    emxArray_creal_T c_v;
    uint32_T d_loop_ub;
    uint32_T k;
    creal_T b[183];
    creal_T d_v[190];
    static const real_T dv0[8] = { 0.16290171402561984, 0.5054728575456503, 0.44610006912318972, -0.019787513117909966, -0.13225358368436993, 0.021808150237390023, 0.023251800535559971, -0.0074934946651300143 };
    creal_T b_b[263];
    creal_T e_v[270];
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'nondecimatedWaveletReconstruction3D:25' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:26' q=length(w1); */
    /* 'nondecimatedWaveletReconstruction3D:27' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    i2 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 270;
    B->size[1] = 190;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i2, (int32_T)sizeof(creal_T));
    loop_ub = 51300 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
        B->data[i2].re = 0.0;
        B->data[i2].im = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:28' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r5, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:29' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletReconstruction3D:30' v=conv(squeeze(A(i, j, :)), w1); */
            i2 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i2, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i2 = 0; i2 <= loop_ub; i2++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i2] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i2];
            }
            squeeze(b_A, v);
            i2 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i2, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i2 = 0; i2 <= loop_ub; i2++) {
                b_v->data[i2] = v->data[i2];
            }
            conv(b_v, v);
            /* 'nondecimatedWaveletReconstruction3D:31' B(i, j, :)=v; */
            i2 = r5->size[0];
            r5->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r5, i2, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (i2 = 0; i2 <= loop_ub; i2++) {
                r5->data[i2] = i2;
            }
            iv0[0] = 1;
            iv0[1] = 1;
            iv0[2] = r5->size[0];
            loop_ub = iv0[2] - 1;
            for (i2 = 0; i2 <= loop_ub; i2++) {
                b_loop_ub = iv0[1] - 1;
                for (i3 = 0; i3 <= b_loop_ub; i3++) {
                    c_loop_ub = iv0[0] - 1;
                    for (i4 = 0; i4 <= c_loop_ub; i4++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv0;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r5->data[i2]] = c_v.data[(i4 + c_v.size[0] * i3) + c_v.size[0] * c_v.size[1] * i2];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r5);
    emxFree_creal_T(&v);
    /* 'nondecimatedWaveletReconstruction3D:35' for i=1:m */
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletReconstruction3D:36' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:37' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (loop_ub = 0; loop_ub < 183; loop_ub++) {
                b[loop_ub] = B->data[(i + B->size[0] * loop_ub) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv0, d_v);
            /* 'nondecimatedWaveletReconstruction3D:38' B(i, :, k)=v; */
            loop_ub = (int32_T)k;
            for (i2 = 0; i2 < 190; i2++) {
                B->data[(i + B->size[0] * i2) + B->size[0] * B->size[1] * (loop_ub - 1)] = d_v[i2];
            }
        }
    }
    /* 'nondecimatedWaveletReconstruction3D:42' for j=1:n+q-1 */
    for (j = 0; j < 190; j++) {
        /* 'nondecimatedWaveletReconstruction3D:43' for k=1:p+q-1 */
        d_loop_ub = (uint32_T)p + 7U;
        for (k = 1U; k <= d_loop_ub; k++) {
            /* 'nondecimatedWaveletReconstruction3D:44' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (loop_ub = 0; loop_ub < 263; loop_ub++) {
                b_b[loop_ub] = B->data[(loop_ub + B->size[0] * j) + B->size[0] * B->size[1] * ((int32_T)k - 1)];
            }
            c_conv(b_b, dv0, e_v);
            /* 'nondecimatedWaveletReconstruction3D:45' B(:, j, k)=v; */
            loop_ub = (int32_T)k;
            for (i2 = 0; i2 < 270; i2++) {
                B->data[(i2 + B->size[0] * j) + B->size[0] * B->size[1] * (loop_ub - 1)] = e_v[i2];
            }
        }
    }
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * 
 */
static void squeeze(const emxArray_creal_T *a, emxArray_creal_T *b)
{
    int32_T k;
    int32_T loop_ub;
    int32_T sqsz[3];
    k = 3;
    while ((k > 2) && (a->size[2] == 1)) {
        k = 2;
    }
    if (k <= 2) {
        k = b->size[0];
        b->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)b, k, (int32_T)sizeof(creal_T));
        loop_ub = a->size[2];
        k = 1;
        while (k <= loop_ub) {
            b->data[0] = a->data[k - 1];
            k = 2;
        }
    } else {
        for (k = 0; k < 3; k++) {
            sqsz[k] = 1;
        }
        if (a->size[2] != 1) {
            sqsz[0] = a->size[2];
        }
        sqsz[1] = 1;
        sqsz[2] = 1;
        k = b->size[0];
        b->size[0] = sqsz[0];
        emxEnsureCapacity((emxArray__common *)b, k, (int32_T)sizeof(creal_T));
        loop_ub = a->size[2];
        for (k = 0; k + 1 <= loop_ub; k++) {
            b->data[k] = a->data[k];
        }
    }
}

/*
 * function B=nondecimatedWaveletReconstruction3D(A)
 */
void nondecimatedWaveletReconstruction3D(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    emxArray_real_T *b_B;
    int32_T p;
    real_T b_p;
    real_T L;
    int32_T loop_ub;
    emxArray_creal_T *b_A;
    int32_T c_B;
    int32_T i0;
    emxArray_creal_T *d_B;
    real_T d0;
    real_T d1;
    emxArray_creal_T *c_A;
    int32_T i1;
    emxArray_creal_T *r0;
    emxArray_creal_T *d_A;
    emxArray_creal_T *r1;
    emxArray_creal_T *e_A;
    emxArray_creal_T *r2;
    emxArray_creal_T *f_A;
    emxArray_creal_T *r3;
    emxArray_creal_T *g_A;
    emxArray_creal_T *r4;
    emxArray_creal_T *h_A;
    emxArray_creal_T *i_A;
    emlrtHeapReferenceStackEnterFcn();
    emxInit_real_T(&b_B, 3, TRUE);
    /* 'nondecimatedWaveletReconstruction3D:2' LR=[2.30377813308855e-001 7.14846570552542e-001 6.30880767929590e-001 -2.79837694169838e-002 -1.87034811718881e-001 3.08413818359870e-002 3.28830116669829e-002 -1.05974017849973e-002]; */
    /* 'nondecimatedWaveletReconstruction3D:3' HR=[-1.05974017849973e-002 -3.28830116669829e-002 3.08413818359870e-002 1.87034811718881e-001 -2.79837694169838e-002 -6.30880767929590e-001 7.14846570552542e-001 -2.30377813308855e-001]; */
    /* 'nondecimatedWaveletReconstruction3D:4' alpha=sum(LR); */
    /* 'nondecimatedWaveletReconstruction3D:4' LR=LR/alpha; */
    /* 'nondecimatedWaveletReconstruction3D:4' HR=HR/alpha; */
    /* 'nondecimatedWaveletReconstruction3D:5' B=waveletReconstruction3D(A, LR, HR); */
    /* 'nondecimatedWaveletReconstruction3D:9' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletReconstruction3D:9' p=p/8; */
    b_p = (real_T)p / 8.0;
    /* 'nondecimatedWaveletReconstruction3D:10' q=length(Lo_R); */
    /* 'nondecimatedWaveletReconstruction3D:10' L=size(A, 3)/8; */
    L = (real_T)A->size[2] / 8.0;
    /* 'nondecimatedWaveletReconstruction3D:11' B=zeros(m+q-1, n+q-1, p+q-1); */
    p = b_B->size[0] * b_B->size[1] * b_B->size[2];
    b_B->size[0] = 270;
    b_B->size[1] = 190;
    b_B->size[2] = (int32_T)emlrtIntegerCheckR2011a((b_p + 8.0) - 1.0, &q_emlrtDCI, &emlrtContextGlobal);
    emxEnsureCapacity((emxArray__common *)b_B, p, (int32_T)sizeof(real_T));
    loop_ub = 51300 * (int32_T)emlrtIntegerCheckR2011a((b_p + 8.0) - 1.0, &r_emlrtDCI, &emlrtContextGlobal) - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_B->data[p] = 0.0;
    }
    /* 'nondecimatedWaveletReconstruction3D:12' B=B+singleBlockRecon(A(:, :, 1:L), Lo_R, Lo_R, Lo_R); */
    if (1.0 > L) {
        p = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(L, &emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&b_A, 3, TRUE);
    c_B = b_A->size[0] * b_A->size[1] * b_A->size[2];
    b_A->size[0] = 263;
    b_A->size[1] = 183;
    b_A->size[2] = p;
    emxEnsureCapacity((emxArray__common *)b_A, c_B, (int32_T)sizeof(creal_T));
    loop_ub = p - 1;
    for (p = 0; p <= loop_ub; p++) {
        for (c_B = 0; c_B < 183; c_B++) {
            for (i0 = 0; i0 < 263; i0++) {
                b_A->data[(i0 + b_A->size[0] * c_B) + b_A->size[0] * b_A->size[1] * p] = A->data[(i0 + A->size[0] * c_B) + A->size[0] * A->size[1] * p];
            }
        }
    }
    emxInit_creal_T(&d_B, 3, TRUE);
    singleBlockRecon(b_A, d_B);
    /* 'nondecimatedWaveletReconstruction3D:13' B=B+singleBlockRecon(A(:, :, L+1:2*L), Hi_R, Lo_R, Lo_R); */
    d0 = L + 1.0;
    d1 = 2.0 * L;
    emxFree_creal_T(&b_A);
    if (d0 > d1) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d0, &b_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d1, &c_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&c_A, 3, TRUE);
    i0 = c_A->size[0] * c_A->size[1] * c_A->size[2];
    c_A->size[0] = 263;
    c_A->size[1] = 183;
    c_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)c_A, i0, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i0 = 0; i0 < 183; i0++) {
            for (i1 = 0; i1 < 263; i1++) {
                c_A->data[(i1 + c_A->size[0] * i0) + c_A->size[0] * c_A->size[1] * c_B] = A->data[(i1 + A->size[0] * i0) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r0, 3, TRUE);
    b_singleBlockRecon(c_A, r0);
    /* 'nondecimatedWaveletReconstruction3D:14' B=B+singleBlockRecon(A(:, :, 2*L+1:3*L), Lo_R, Hi_R, Lo_R); */
    d0 = 2.0 * L + 1.0;
    d1 = 3.0 * L;
    emxFree_creal_T(&c_A);
    if (d0 > d1) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d0, &d_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d1, &e_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&d_A, 3, TRUE);
    i0 = d_A->size[0] * d_A->size[1] * d_A->size[2];
    d_A->size[0] = 263;
    d_A->size[1] = 183;
    d_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)d_A, i0, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i0 = 0; i0 < 183; i0++) {
            for (i1 = 0; i1 < 263; i1++) {
                d_A->data[(i1 + d_A->size[0] * i0) + d_A->size[0] * d_A->size[1] * c_B] = A->data[(i1 + A->size[0] * i0) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r1, 3, TRUE);
    c_singleBlockRecon(d_A, r1);
    /* 'nondecimatedWaveletReconstruction3D:15' B=B+singleBlockRecon(A(:, :, 3*L+1:4*L), Hi_R, Hi_R, Lo_R); */
    d0 = 3.0 * L + 1.0;
    d1 = 4.0 * L;
    emxFree_creal_T(&d_A);
    if (d0 > d1) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d0, &f_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d1, &g_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&e_A, 3, TRUE);
    i0 = e_A->size[0] * e_A->size[1] * e_A->size[2];
    e_A->size[0] = 263;
    e_A->size[1] = 183;
    e_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)e_A, i0, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i0 = 0; i0 < 183; i0++) {
            for (i1 = 0; i1 < 263; i1++) {
                e_A->data[(i1 + e_A->size[0] * i0) + e_A->size[0] * e_A->size[1] * c_B] = A->data[(i1 + A->size[0] * i0) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r2, 3, TRUE);
    d_singleBlockRecon(e_A, r2);
    /* 'nondecimatedWaveletReconstruction3D:16' B=B+singleBlockRecon(A(:, :, 4*L+1:5*L), Lo_R, Lo_R, Hi_R); */
    d0 = 4.0 * L + 1.0;
    d1 = 5.0 * L;
    emxFree_creal_T(&e_A);
    if (d0 > d1) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d0, &h_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d1, &i_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&f_A, 3, TRUE);
    i0 = f_A->size[0] * f_A->size[1] * f_A->size[2];
    f_A->size[0] = 263;
    f_A->size[1] = 183;
    f_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)f_A, i0, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i0 = 0; i0 < 183; i0++) {
            for (i1 = 0; i1 < 263; i1++) {
                f_A->data[(i1 + f_A->size[0] * i0) + f_A->size[0] * f_A->size[1] * c_B] = A->data[(i1 + A->size[0] * i0) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r3, 3, TRUE);
    e_singleBlockRecon(f_A, r3);
    /* 'nondecimatedWaveletReconstruction3D:17' B=B+singleBlockRecon(A(:, :, 5*L+1:6*L), Hi_R, Lo_R, Hi_R); */
    d0 = 5.0 * L + 1.0;
    d1 = 6.0 * L;
    emxFree_creal_T(&f_A);
    if (d0 > d1) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d0, &j_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d1, &k_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&g_A, 3, TRUE);
    i0 = g_A->size[0] * g_A->size[1] * g_A->size[2];
    g_A->size[0] = 263;
    g_A->size[1] = 183;
    g_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)g_A, i0, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i0 = 0; i0 < 183; i0++) {
            for (i1 = 0; i1 < 263; i1++) {
                g_A->data[(i1 + g_A->size[0] * i0) + g_A->size[0] * g_A->size[1] * c_B] = A->data[(i1 + A->size[0] * i0) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r4, 3, TRUE);
    f_singleBlockRecon(g_A, r4);
    p = d_B->size[0] * d_B->size[1] * d_B->size[2];
    d_B->size[0] = 270;
    d_B->size[1] = 190;
    d_B->size[2] = b_B->size[2];
    emxEnsureCapacity((emxArray__common *)d_B, p, (int32_T)sizeof(creal_T));
    emxFree_creal_T(&g_A);
    loop_ub = b_B->size[0] * b_B->size[1] * b_B->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        d_B->data[p].re = (((((b_B->data[p] + d_B->data[p].re) + r0->data[p].re) + r1->data[p].re) + r2->data[p].re) + r3->data[p].re) + r4->data[p].re;
        d_B->data[p].im = ((((d_B->data[p].im + r0->data[p].im) + r1->data[p].im) + r2->data[p].im) + r3->data[p].im) + r4->data[p].im;
    }
    emxFree_creal_T(&r4);
    emxFree_creal_T(&r3);
    emxFree_creal_T(&r2);
    emxFree_creal_T(&r1);
    emxFree_real_T(&b_B);
    /* 'nondecimatedWaveletReconstruction3D:18' B=B+singleBlockRecon(A(:, :, 6*L+1:7*L), Lo_R, Hi_R, Hi_R); */
    d0 = 6.0 * L + 1.0;
    d1 = 7.0 * L;
    if (d0 > d1) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d0, &l_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d1, &m_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&h_A, 3, TRUE);
    i0 = h_A->size[0] * h_A->size[1] * h_A->size[2];
    h_A->size[0] = 263;
    h_A->size[1] = 183;
    h_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)h_A, i0, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i0 = 0; i0 < 183; i0++) {
            for (i1 = 0; i1 < 263; i1++) {
                h_A->data[(i1 + h_A->size[0] * i0) + h_A->size[0] * h_A->size[1] * c_B] = A->data[(i1 + A->size[0] * i0) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    g_singleBlockRecon(h_A, r0);
    p = d_B->size[0] * d_B->size[1] * d_B->size[2];
    d_B->size[0] = 270;
    d_B->size[1] = 190;
    emxEnsureCapacity((emxArray__common *)d_B, p, (int32_T)sizeof(creal_T));
    p = d_B->size[0];
    loop_ub = d_B->size[1];
    c_B = d_B->size[2];
    emxFree_creal_T(&h_A);
    loop_ub = p * loop_ub * c_B - 1;
    for (p = 0; p <= loop_ub; p++) {
        d_B->data[p].re += r0->data[p].re;
        d_B->data[p].im += r0->data[p].im;
    }
    /* 'nondecimatedWaveletReconstruction3D:19' B=B+singleBlockRecon(A(:, :, 7*L+1:8*L), Hi_R, Hi_R, Hi_R); */
    d0 = 7.0 * L + 1.0;
    d1 = 8.0 * L;
    if (d0 > d1) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d0, &n_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d1, &o_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&i_A, 3, TRUE);
    i0 = i_A->size[0] * i_A->size[1] * i_A->size[2];
    i_A->size[0] = 263;
    i_A->size[1] = 183;
    i_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)i_A, i0, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i0 = 0; i0 < 183; i0++) {
            for (i1 = 0; i1 < 263; i1++) {
                i_A->data[(i1 + i_A->size[0] * i0) + i_A->size[0] * i_A->size[1] * c_B] = A->data[(i1 + A->size[0] * i0) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    h_singleBlockRecon(i_A, r0);
    p = d_B->size[0] * d_B->size[1] * d_B->size[2];
    d_B->size[0] = 270;
    d_B->size[1] = 190;
    emxEnsureCapacity((emxArray__common *)d_B, p, (int32_T)sizeof(creal_T));
    p = d_B->size[0];
    loop_ub = d_B->size[1];
    c_B = d_B->size[2];
    emxFree_creal_T(&i_A);
    loop_ub = p * loop_ub * c_B - 1;
    for (p = 0; p <= loop_ub; p++) {
        d_B->data[p].re += r0->data[p].re;
        d_B->data[p].im += r0->data[p].im;
    }
    emxFree_creal_T(&r0);
    /* 'nondecimatedWaveletReconstruction3D:20' B=B(q:m, q:n, q:p); */
    if (8.0 > b_p) {
        p = 0;
        c_B = 0;
    } else {
        p = 7;
        c_B = (int32_T)emlrtIntegerCheckR2011a(b_p, &p_emlrtDCI, &emlrtContextGlobal);
    }
    i0 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 256;
    B->size[1] = 176;
    B->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)B, i0, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i0 = 0; i0 < 176; i0++) {
            for (i1 = 0; i1 < 256; i1++) {
                B->data[(i1 + B->size[0] * i0) + B->size[0] * B->size[1] * c_B] = d_B->data[((i1 + d_B->size[0] * (7 + i0)) + d_B->size[0] * d_B->size[1] * (p + c_B)) + 7];
            }
        }
    }
    emxFree_creal_T(&d_B);
    emlrtHeapReferenceStackLeaveFcn();
}

void nondecimatedWaveletReconstruction3D_api(const mxArray * const prhs[1], const mxArray *plhs[1])
{
    emxArray_creal_T *A;
    emxArray_creal_T *B;
    emlrtHeapReferenceStackEnterFcn();
    emxInit_creal_T(&A, 3, TRUE);
    emxInit_creal_T(&B, 3, TRUE);
    /* Marshall function inputs */
    emlrt_marshallIn(emlrtAliasP(prhs[0]), "A", A);
    /* Invoke the target function */
    nondecimatedWaveletReconstruction3D(A, B);
    /* Marshall function outputs */
    plhs[0] = emlrt_marshallOut(B);
    emxFree_creal_T(&B);
    emxFree_creal_T(&A);
    emlrtHeapReferenceStackLeaveFcn();
}

void nondecimatedWaveletReconstruction3D_atexit(void)
{
    emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void nondecimatedWaveletReconstruction3D_initialize(emlrtContext *context)
{
    emlrtEnterRtStack(&emlrtContextGlobal);
    emlrtFirstTime(context);
}

void nondecimatedWaveletReconstruction3D_terminate(void)
{
    emlrtLeaveRtStack(&emlrtContextGlobal);
}
/* End of code generation (nondecimatedWaveletReconstruction3D.c) */
