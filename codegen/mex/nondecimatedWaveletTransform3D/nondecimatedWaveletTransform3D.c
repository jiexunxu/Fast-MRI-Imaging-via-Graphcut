/*
 * nondecimatedWaveletTransform3D.c
 *
 * Code generation for function 'nondecimatedWaveletTransform3D'
 *
 * C source code generated on: Thu Sep 06 12:09:03 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nondecimatedWaveletTransform3D.h"

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

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_conv(const creal_T A[176], const real_T B[8], creal_T C[183]);
static void b_directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B);
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_creal_T *y);
static void b_emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void b_squeeze(const emxArray_creal_T *a, emxArray_creal_T *b);
static void c_directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B);
static void c_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_creal_T *ret);
static void conv(const creal_T A[256], const real_T B[8], creal_T C[263]);
static void d_directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B);
static void directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B);
static void e_directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B);
static void emlrt_marshallIn(const mxArray *A, const char_T *identifier, emxArray_creal_T *y);
static const mxArray *emlrt_marshallOut(emxArray_creal_T *u);
static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize);
static void emxFree_creal_T(emxArray_creal_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void squeeze(const creal_T a[256], creal_T b[256]);

/* Function Definitions */

/*
 * 
 */
static void b_conv(const creal_T A[176], const real_T B[8], creal_T C[183])
{
    int32_T jC;
    int32_T jA1;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    for (jC = 0; jC < 183; jC++) {
        if (8 < jC + 2) {
            jA1 = jC;
        } else {
            jA1 = 7;
        }
        if (176 < jC + 1) {
            jA2 = 176;
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

/*
 * function B=directionalConv(A, w, dir)
 */
static void b_directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i2;
    int32_T loop_ub;
    uint32_T k;
    int32_T b_k;
    creal_T b[176];
    creal_T v[183];
    static const real_T dv2[8] = { -0.0074934946651300152, 0.023251800535559974, 0.021808150237390026, -0.13225358368436993, -0.01978751311790997, 0.44610006912318978, 0.5054728575456503, 0.16290171402561987 };
    /* 'nondecimatedWaveletTransform3D:31' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletTransform3D:32' q=length(w); */
    /* 'nondecimatedWaveletTransform3D:33' if dir==1 */
    /* 'nondecimatedWaveletTransform3D:41' elseif dir==2 */
    /* 'nondecimatedWaveletTransform3D:42' B=complex(zeros(m, n+q-1, p), 0); */
    i2 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 263;
    B->size[1] = 183;
    B->size[2] = p;
    emxEnsureCapacity((emxArray__common *)B, i2, (int32_T)sizeof(creal_T));
    loop_ub = 48129 * p - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
        B->data[i2].re = 0.0;
        B->data[i2].im = 0.0;
    }
    /* 'nondecimatedWaveletTransform3D:43' for i=1:m */
    for (loop_ub = 0; loop_ub < 263; loop_ub++) {
        /* 'nondecimatedWaveletTransform3D:44' for k=1:p */
        for (k = 1U; k <= (uint32_T)p; k++) {
            /* 'nondecimatedWaveletTransform3D:45' v=conv(squeeze(A(i, :, k)), w); */
            for (b_k = 0; b_k < 176; b_k++) {
                b[b_k] = A->data[(loop_ub + A->size[0] * b_k) + A->size[0] * A->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv2, v);
            /* 'nondecimatedWaveletTransform3D:46' B(i, :, k)=v; */
            b_k = (int32_T)k;
            for (i2 = 0; i2 < 183; i2++) {
                B->data[(loop_ub + B->size[0] * i2) + B->size[0] * B->size[1] * (b_k - 1)] = v[i2];
            }
        }
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
 * 
 */
static void b_squeeze(const emxArray_creal_T *a, emxArray_creal_T *b)
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
 * function B=directionalConv(A, w, dir)
 */
static void c_directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i3;
    int32_T loop_ub;
    uint32_T k;
    int32_T b_k;
    creal_T b[176];
    creal_T v[183];
    static const real_T dv3[8] = { -0.16290171402561987, 0.5054728575456503, -0.44610006912318978, -0.01978751311790997, 0.13225358368436993, 0.021808150237390026, -0.023251800535559974, -0.0074934946651300152 };
    /* 'nondecimatedWaveletTransform3D:31' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletTransform3D:32' q=length(w); */
    /* 'nondecimatedWaveletTransform3D:33' if dir==1 */
    /* 'nondecimatedWaveletTransform3D:41' elseif dir==2 */
    /* 'nondecimatedWaveletTransform3D:42' B=complex(zeros(m, n+q-1, p), 0); */
    i3 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 263;
    B->size[1] = 183;
    B->size[2] = p;
    emxEnsureCapacity((emxArray__common *)B, i3, (int32_T)sizeof(creal_T));
    loop_ub = 48129 * p - 1;
    for (i3 = 0; i3 <= loop_ub; i3++) {
        B->data[i3].re = 0.0;
        B->data[i3].im = 0.0;
    }
    /* 'nondecimatedWaveletTransform3D:43' for i=1:m */
    for (loop_ub = 0; loop_ub < 263; loop_ub++) {
        /* 'nondecimatedWaveletTransform3D:44' for k=1:p */
        for (k = 1U; k <= (uint32_T)p; k++) {
            /* 'nondecimatedWaveletTransform3D:45' v=conv(squeeze(A(i, :, k)), w); */
            for (b_k = 0; b_k < 176; b_k++) {
                b[b_k] = A->data[(loop_ub + A->size[0] * b_k) + A->size[0] * A->size[1] * ((int32_T)k - 1)];
            }
            b_conv(b, dv3, v);
            /* 'nondecimatedWaveletTransform3D:46' B(i, :, k)=v; */
            b_k = (int32_T)k;
            for (i3 = 0; i3 < 183; i3++) {
                B->data[(loop_ub + B->size[0] * i3) + B->size[0] * B->size[1] * (b_k - 1)] = v[i3];
            }
        }
    }
}

static void c_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_creal_T *ret)
{
    int32_T i;
    static const int16_T iv2[3] = { 256, 176, 256 };
    int32_T iv3[3];
    static const boolean_T bv0[3] = { FALSE, FALSE, TRUE };
    boolean_T bv1[3];
    for (i = 0; i < 3; i++) {
        iv3[i] = iv2[i];
        bv1[i] = bv0[i];
    }
    emlrtCheckVsBuiltInR2011a(msgId, src, "double", TRUE, 3U, iv3, bv1, ret->size);
    i = ret->size[0] * ret->size[1] * ret->size[2];
    ret->size[0] = 256;
    ret->size[1] = 176;
    emxEnsureCapacity((emxArray__common *)ret, i, (int32_T)sizeof(creal_T));
    emlrtImportArrayR2008b(src, ret->data, 8);
    emlrtDestroyArray(&src);
}

/*
 * 
 */
static void conv(const creal_T A[256], const real_T B[8], creal_T C[263])
{
    int32_T jC;
    int32_T jA1;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    for (jC = 0; jC < 263; jC++) {
        if (8 < jC + 2) {
            jA1 = jC;
        } else {
            jA1 = 7;
        }
        if (256 < jC + 1) {
            jA2 = 256;
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

/*
 * function B=directionalConv(A, w, dir)
 */
static void d_directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i4;
    emxArray_creal_T *v;
    emxArray_int32_T *r1;
    emxArray_creal_T *b_A;
    emxArray_creal_T *c_A;
    int32_T i;
    int32_T j;
    int32_T nA;
    int32_T nApnB;
    int32_T jA1;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    static const real_T dv4[8] = { -0.0074934946651300152, 0.023251800535559974, 0.021808150237390026, -0.13225358368436993, -0.01978751311790997, 0.44610006912318978, 0.5054728575456503, 0.16290171402561987 };
    int32_T iv0[3];
    emxArray_creal_T b_v;
    emlrtHeapReferenceStackEnterFcn();
    /* 'nondecimatedWaveletTransform3D:31' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletTransform3D:32' q=length(w); */
    /* 'nondecimatedWaveletTransform3D:33' if dir==1 */
    /* 'nondecimatedWaveletTransform3D:49' elseif dir==3 */
    /* 'nondecimatedWaveletTransform3D:50' B=complex(zeros(m, n, p+q-1), 0); */
    i4 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 263;
    B->size[1] = 183;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i4, (int32_T)sizeof(creal_T));
    p = 48129 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i4 = 0; i4 <= p; i4++) {
        B->data[i4].re = 0.0;
        B->data[i4].im = 0.0;
    }
    /* 'nondecimatedWaveletTransform3D:51' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r1, 1, TRUE);
    b_emxInit_creal_T(&b_A, 1, TRUE);
    emxInit_creal_T(&c_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletTransform3D:52' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletTransform3D:53' v=conv(squeeze(A(i, j, :)), w); */
            i4 = c_A->size[0] * c_A->size[1] * c_A->size[2];
            c_A->size[0] = 1;
            c_A->size[1] = 1;
            c_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)c_A, i4, (int32_T)sizeof(creal_T));
            p = A->size[2] - 1;
            for (i4 = 0; i4 <= p; i4++) {
                c_A->data[c_A->size[0] * c_A->size[1] * i4] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i4];
            }
            b_squeeze(c_A, b_A);
            nA = b_A->size[0];
            nApnB = nA + 7;
            if (nA == 0) {
                nApnB++;
            }
            i4 = v->size[0];
            v->size[0] = nApnB;
            emxEnsureCapacity((emxArray__common *)v, i4, (int32_T)sizeof(creal_T));
            if ((b_A->size[0] == 0) || (nApnB == 0)) {
                p = v->size[0];
                i4 = v->size[0];
                v->size[0] = p;
                emxEnsureCapacity((emxArray__common *)v, i4, (int32_T)sizeof(creal_T));
                p--;
                for (i4 = 0; i4 <= p; i4++) {
                    v->data[i4].re = 0.0;
                    v->data[i4].im = 0.0;
                }
            } else {
                for (p = 1; p <= nApnB; p++) {
                    if (8 < p + 1) {
                        jA1 = p - 7;
                    } else {
                        jA1 = 1;
                    }
                    if (nA < p) {
                        jA2 = nA;
                    } else {
                        jA2 = p;
                    }
                    s_re = 0.0;
                    s_im = 0.0;
                    while (jA1 <= jA2) {
                        s_re += b_A->data[jA1 - 1].re * dv4[p - jA1];
                        s_im += b_A->data[jA1 - 1].im * dv4[p - jA1];
                        jA1++;
                    }
                    v->data[p - 1].re = s_re;
                    v->data[p - 1].im = s_im;
                }
            }
            /* 'nondecimatedWaveletTransform3D:54' B(i, j, :)=v; */
            i4 = r1->size[0];
            r1->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r1, i4, (int32_T)sizeof(int32_T));
            p = B->size[2] - 1;
            for (i4 = 0; i4 <= p; i4++) {
                r1->data[i4] = i4;
            }
            iv0[0] = 1;
            iv0[1] = 1;
            iv0[2] = r1->size[0];
            p = iv0[2] - 1;
            for (i4 = 0; i4 <= p; i4++) {
                jA1 = iv0[1] - 1;
                for (jA2 = 0; jA2 <= jA1; jA2++) {
                    nA = iv0[0] - 1;
                    for (nApnB = 0; nApnB <= nA; nApnB++) {
                        b_v = *v;
                        b_v.size = (int32_T *)&iv0;
                        b_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r1->data[i4]] = b_v.data[(nApnB + b_v.size[0] * jA2) + b_v.size[0] * b_v.size[1] * i4];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&c_A);
    emxFree_creal_T(&b_A);
    emxFree_int32_T(&r1);
    emxFree_creal_T(&v);
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * function B=directionalConv(A, w, dir)
 */
static void directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i1;
    int32_T loop_ub;
    uint32_T k;
    creal_T dcv1[256];
    creal_T v[263];
    static const real_T dv1[8] = { -0.0074934946651300152, 0.023251800535559974, 0.021808150237390026, -0.13225358368436993, -0.01978751311790997, 0.44610006912318978, 0.5054728575456503, 0.16290171402561987 };
    int32_T b_k;
    /* 'nondecimatedWaveletTransform3D:31' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletTransform3D:32' q=length(w); */
    /* 'nondecimatedWaveletTransform3D:33' if dir==1 */
    /* 'nondecimatedWaveletTransform3D:34' B=complex(zeros(m+q-1, n, p), 0); */
    i1 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 263;
    B->size[1] = 176;
    B->size[2] = p;
    emxEnsureCapacity((emxArray__common *)B, i1, (int32_T)sizeof(creal_T));
    loop_ub = 46288 * p - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
        B->data[i1].re = 0.0;
        B->data[i1].im = 0.0;
    }
    /* 'nondecimatedWaveletTransform3D:35' for j=1:n */
    for (loop_ub = 0; loop_ub < 176; loop_ub++) {
        /* 'nondecimatedWaveletTransform3D:36' for k=1:p */
        for (k = 1U; k <= (uint32_T)p; k++) {
            /* 'nondecimatedWaveletTransform3D:37' v=conv(squeeze(A(:, j, k)), w); */
            squeeze(*(creal_T (*)[256])&A->data[A->size[0] * loop_ub + A->size[0] * A->size[1] * ((int32_T)k - 1)], dcv1);
            conv(dcv1, dv1, v);
            /* 'nondecimatedWaveletTransform3D:38' B(:, j, k)=v; */
            b_k = (int32_T)k;
            for (i1 = 0; i1 < 263; i1++) {
                B->data[(i1 + B->size[0] * loop_ub) + B->size[0] * B->size[1] * (b_k - 1)] = v[i1];
            }
        }
    }
}

/*
 * function B=directionalConv(A, w, dir)
 */
static void e_directionalConv(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T i5;
    emxArray_creal_T *v;
    emxArray_int32_T *r2;
    emxArray_creal_T *b_A;
    emxArray_creal_T *c_A;
    int32_T i;
    int32_T j;
    int32_T nA;
    int32_T nApnB;
    int32_T jA1;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    static const real_T dv5[8] = { -0.16290171402561987, 0.5054728575456503, -0.44610006912318978, -0.01978751311790997, 0.13225358368436993, 0.021808150237390026, -0.023251800535559974, -0.0074934946651300152 };
    int32_T iv1[3];
    emxArray_creal_T b_v;
    emlrtHeapReferenceStackEnterFcn();
    /* 'nondecimatedWaveletTransform3D:31' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'nondecimatedWaveletTransform3D:32' q=length(w); */
    /* 'nondecimatedWaveletTransform3D:33' if dir==1 */
    /* 'nondecimatedWaveletTransform3D:49' elseif dir==3 */
    /* 'nondecimatedWaveletTransform3D:50' B=complex(zeros(m, n, p+q-1), 0); */
    i5 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 263;
    B->size[1] = 183;
    B->size[2] = (int32_T)(((real_T)p + 8.0) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i5, (int32_T)sizeof(creal_T));
    p = 48129 * (int32_T)(((real_T)p + 8.0) - 1.0) - 1;
    for (i5 = 0; i5 <= p; i5++) {
        B->data[i5].re = 0.0;
        B->data[i5].im = 0.0;
    }
    /* 'nondecimatedWaveletTransform3D:51' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r2, 1, TRUE);
    b_emxInit_creal_T(&b_A, 1, TRUE);
    emxInit_creal_T(&c_A, 3, TRUE);
    for (i = 0; i < 263; i++) {
        /* 'nondecimatedWaveletTransform3D:52' for j=1:n */
        for (j = 0; j < 183; j++) {
            /* 'nondecimatedWaveletTransform3D:53' v=conv(squeeze(A(i, j, :)), w); */
            i5 = c_A->size[0] * c_A->size[1] * c_A->size[2];
            c_A->size[0] = 1;
            c_A->size[1] = 1;
            c_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)c_A, i5, (int32_T)sizeof(creal_T));
            p = A->size[2] - 1;
            for (i5 = 0; i5 <= p; i5++) {
                c_A->data[c_A->size[0] * c_A->size[1] * i5] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * i5];
            }
            b_squeeze(c_A, b_A);
            nA = b_A->size[0];
            nApnB = nA + 7;
            if (nA == 0) {
                nApnB++;
            }
            i5 = v->size[0];
            v->size[0] = nApnB;
            emxEnsureCapacity((emxArray__common *)v, i5, (int32_T)sizeof(creal_T));
            if ((b_A->size[0] == 0) || (nApnB == 0)) {
                p = v->size[0];
                i5 = v->size[0];
                v->size[0] = p;
                emxEnsureCapacity((emxArray__common *)v, i5, (int32_T)sizeof(creal_T));
                p--;
                for (i5 = 0; i5 <= p; i5++) {
                    v->data[i5].re = 0.0;
                    v->data[i5].im = 0.0;
                }
            } else {
                for (p = 1; p <= nApnB; p++) {
                    if (8 < p + 1) {
                        jA1 = p - 7;
                    } else {
                        jA1 = 1;
                    }
                    if (nA < p) {
                        jA2 = nA;
                    } else {
                        jA2 = p;
                    }
                    s_re = 0.0;
                    s_im = 0.0;
                    while (jA1 <= jA2) {
                        s_re += b_A->data[jA1 - 1].re * dv5[p - jA1];
                        s_im += b_A->data[jA1 - 1].im * dv5[p - jA1];
                        jA1++;
                    }
                    v->data[p - 1].re = s_re;
                    v->data[p - 1].im = s_im;
                }
            }
            /* 'nondecimatedWaveletTransform3D:54' B(i, j, :)=v; */
            i5 = r2->size[0];
            r2->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r2, i5, (int32_T)sizeof(int32_T));
            p = B->size[2] - 1;
            for (i5 = 0; i5 <= p; i5++) {
                r2->data[i5] = i5;
            }
            iv1[0] = 1;
            iv1[1] = 1;
            iv1[2] = r2->size[0];
            p = iv1[2] - 1;
            for (i5 = 0; i5 <= p; i5++) {
                jA1 = iv1[1] - 1;
                for (jA2 = 0; jA2 <= jA1; jA2++) {
                    nA = iv1[0] - 1;
                    for (nApnB = 0; nApnB <= nA; nApnB++) {
                        b_v = *v;
                        b_v.size = (int32_T *)&iv1;
                        b_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r2->data[i5]] = b_v.data[(nApnB + b_v.size[0] * jA2) + b_v.size[0] * b_v.size[1] * i5];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&c_A);
    emxFree_creal_T(&b_A);
    emxFree_int32_T(&r2);
    emxFree_creal_T(&v);
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

/*
 * 
 */
static void squeeze(const creal_T a[256], creal_T b[256])
{
    memcpy((void *)&b[0], (void *)&a[0], sizeof(creal_T) << 8);
}

/*
 * function B=nondecimatedWaveletTransform3D(A)
 */
void nondecimatedWaveletTransform3D(const emxArray_creal_T *A, emxArray_creal_T *B)
{
    int32_T p;
    int32_T y;
    int32_T i0;
    int32_T loop_ub;
    emxArray_creal_T *BL;
    emxArray_creal_T *BH;
    int32_T b_p;
    int32_T k;
    creal_T dcv0[256];
    creal_T v[263];
    static const real_T dv0[8] = { -0.16290171402561987, 0.5054728575456503, -0.44610006912318978, -0.01978751311790997, 0.13225358368436993, 0.021808150237390026, -0.023251800535559974, -0.0074934946651300152 };
    emxArray_creal_T *BLL;
    emxArray_creal_T *BLH;
    emxArray_creal_T *BHL;
    emxArray_creal_T *BHH;
    emxArray_creal_T *r0;
    emlrtHeapReferenceStackEnterFcn();
    /* 'nondecimatedWaveletTransform3D:2' LD=[-1.05974017849973e-002 3.28830116669829e-002 3.08413818359870e-002 -1.87034811718881e-001 -2.79837694169838e-002 6.30880767929590e-001 7.14846570552542e-001 2.30377813308855e-001]; */
    /* 'nondecimatedWaveletTransform3D:3' HD=[-2.30377813308855e-001 7.14846570552542e-001 -6.30880767929590e-001 -2.79837694169838e-002 1.87034811718881e-001 3.08413818359870e-002 -3.28830116669829e-002 -1.05974017849973e-002]; */
    /* 'nondecimatedWaveletTransform3D:4' alpha=sum(LD); */
    /* 'nondecimatedWaveletTransform3D:4' LD=LD/alpha; */
    /* 'nondecimatedWaveletTransform3D:4' HD=HD/alpha; */
    /* 'nondecimatedWaveletTransform3D:5' B=waveletDecomposition3D(A, LD, HD); */
    /*  LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH */
    /* 'nondecimatedWaveletTransform3D:10' [m, n, p]=size(A); */
    p = A->size[2] + 7;
    /* 'nondecimatedWaveletTransform3D:11' q=length(Lo_D); */
    /* 'nondecimatedWaveletTransform3D:12' B=complex(zeros(m+q-1, n+q-1, (p+q-1)*8), 0); */
    y = p << 3;
    i0 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = 263;
    B->size[1] = 183;
    B->size[2] = y;
    emxEnsureCapacity((emxArray__common *)B, i0, (int32_T)sizeof(creal_T));
    loop_ub = 48129 * y - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        B->data[i0].re = 0.0;
        B->data[i0].im = 0.0;
    }
    emxInit_creal_T(&BL, 3, TRUE);
    emxInit_creal_T(&BH, 3, TRUE);
    /* 'nondecimatedWaveletTransform3D:13' BL=directionalConv(A, Lo_D, 1); */
    directionalConv(A, BL);
    /* 'nondecimatedWaveletTransform3D:14' BH=directionalConv(A, Hi_D, 1); */
    /* 'nondecimatedWaveletTransform3D:31' [m, n, p]=size(A); */
    b_p = A->size[2];
    /* 'nondecimatedWaveletTransform3D:32' q=length(w); */
    /* 'nondecimatedWaveletTransform3D:33' if dir==1 */
    /* 'nondecimatedWaveletTransform3D:34' B=complex(zeros(m+q-1, n, p), 0); */
    i0 = BH->size[0] * BH->size[1] * BH->size[2];
    BH->size[0] = 263;
    BH->size[1] = 176;
    BH->size[2] = b_p;
    emxEnsureCapacity((emxArray__common *)BH, i0, (int32_T)sizeof(creal_T));
    loop_ub = 46288 * b_p - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        BH->data[i0].re = 0.0;
        BH->data[i0].im = 0.0;
    }
    /* 'nondecimatedWaveletTransform3D:35' for j=1:n */
    for (y = 0; y < 176; y++) {
        /* 'nondecimatedWaveletTransform3D:36' for k=1:p */
        for (k = 0; k + 1 <= b_p; k++) {
            /* 'nondecimatedWaveletTransform3D:37' v=conv(squeeze(A(:, j, k)), w); */
            squeeze(*(creal_T (*)[256])&A->data[A->size[0] * y + A->size[0] * A->size[1] * k], dcv0);
            conv(dcv0, dv0, v);
            /* 'nondecimatedWaveletTransform3D:38' B(:, j, k)=v; */
            for (i0 = 0; i0 < 263; i0++) {
                BH->data[(i0 + BH->size[0] * y) + BH->size[0] * BH->size[1] * k] = v[i0];
            }
        }
    }
    emxInit_creal_T(&BLL, 3, TRUE);
    emxInit_creal_T(&BLH, 3, TRUE);
    emxInit_creal_T(&BHL, 3, TRUE);
    emxInit_creal_T(&BHH, 3, TRUE);
    emxInit_creal_T(&r0, 3, TRUE);
    /* 'nondecimatedWaveletTransform3D:15' BLL=directionalConv(BL, Lo_D, 2); */
    b_directionalConv(BL, BLL);
    /* 'nondecimatedWaveletTransform3D:16' BLH=directionalConv(BL, Hi_D, 2); */
    c_directionalConv(BL, BLH);
    /* 'nondecimatedWaveletTransform3D:17' BHL=directionalConv(BH, Lo_D, 2); */
    b_directionalConv(BH, BHL);
    /* 'nondecimatedWaveletTransform3D:18' BHH=directionalConv(BH, Hi_D, 2); */
    c_directionalConv(BH, BHH);
    /* 'nondecimatedWaveletTransform3D:19' L=p+q-1; */
    /* 'nondecimatedWaveletTransform3D:20' B(:, :, 1:L)=directionalConv(BLL, Lo_D, 3); */
    d_directionalConv(BLL, r0);
    emxFree_creal_T(&BH);
    emxFree_creal_T(&BL);
    loop_ub = r0->size[2] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        for (y = 0; y < 183; y++) {
            for (k = 0; k < 263; k++) {
                B->data[(k + B->size[0] * y) + B->size[0] * B->size[1] * i0] = r0->data[(k + r0->size[0] * y) + r0->size[0] * r0->size[1] * i0];
            }
        }
    }
    /* 'nondecimatedWaveletTransform3D:21' B(:, :, L+1:2*L)=directionalConv(BLL, Hi_D, 3); */
    i0 = p + 1;
    if (i0 > (p << 1)) {
        i0 = 0;
    } else {
        i0--;
    }
    e_directionalConv(BLL, r0);
    emxFree_creal_T(&BLL);
    loop_ub = r0->size[2] - 1;
    for (y = 0; y <= loop_ub; y++) {
        for (k = 0; k < 183; k++) {
            for (b_p = 0; b_p < 263; b_p++) {
                B->data[(b_p + B->size[0] * k) + B->size[0] * B->size[1] * (i0 + y)] = r0->data[(b_p + r0->size[0] * k) + r0->size[0] * r0->size[1] * y];
            }
        }
    }
    /* 'nondecimatedWaveletTransform3D:22' B(:, :, 2*L+1:3*L)=directionalConv(BLH, Lo_D, 3); */
    i0 = (p << 1) + 1;
    if (i0 > 3 * p) {
        i0 = 0;
    } else {
        i0--;
    }
    d_directionalConv(BLH, r0);
    loop_ub = r0->size[2] - 1;
    for (y = 0; y <= loop_ub; y++) {
        for (k = 0; k < 183; k++) {
            for (b_p = 0; b_p < 263; b_p++) {
                B->data[(b_p + B->size[0] * k) + B->size[0] * B->size[1] * (i0 + y)] = r0->data[(b_p + r0->size[0] * k) + r0->size[0] * r0->size[1] * y];
            }
        }
    }
    /* 'nondecimatedWaveletTransform3D:23' B(:, :, 3*L+1:4*L)=directionalConv(BLH, Hi_D, 3); */
    i0 = 3 * p + 1;
    if (i0 > (p << 2)) {
        i0 = 0;
    } else {
        i0--;
    }
    e_directionalConv(BLH, r0);
    emxFree_creal_T(&BLH);
    loop_ub = r0->size[2] - 1;
    for (y = 0; y <= loop_ub; y++) {
        for (k = 0; k < 183; k++) {
            for (b_p = 0; b_p < 263; b_p++) {
                B->data[(b_p + B->size[0] * k) + B->size[0] * B->size[1] * (i0 + y)] = r0->data[(b_p + r0->size[0] * k) + r0->size[0] * r0->size[1] * y];
            }
        }
    }
    /* 'nondecimatedWaveletTransform3D:24' B(:, :, 4*L+1:5*L)=directionalConv(BHL, Lo_D, 3); */
    i0 = (p << 2) + 1;
    if (i0 > 5 * p) {
        i0 = 0;
    } else {
        i0--;
    }
    d_directionalConv(BHL, r0);
    loop_ub = r0->size[2] - 1;
    for (y = 0; y <= loop_ub; y++) {
        for (k = 0; k < 183; k++) {
            for (b_p = 0; b_p < 263; b_p++) {
                B->data[(b_p + B->size[0] * k) + B->size[0] * B->size[1] * (i0 + y)] = r0->data[(b_p + r0->size[0] * k) + r0->size[0] * r0->size[1] * y];
            }
        }
    }
    /* 'nondecimatedWaveletTransform3D:25' B(:, :, 5*L+1:6*L)=directionalConv(BHL, Hi_D, 3); */
    i0 = 5 * p + 1;
    if (i0 > 6 * p) {
        i0 = 0;
    } else {
        i0--;
    }
    e_directionalConv(BHL, r0);
    emxFree_creal_T(&BHL);
    loop_ub = r0->size[2] - 1;
    for (y = 0; y <= loop_ub; y++) {
        for (k = 0; k < 183; k++) {
            for (b_p = 0; b_p < 263; b_p++) {
                B->data[(b_p + B->size[0] * k) + B->size[0] * B->size[1] * (i0 + y)] = r0->data[(b_p + r0->size[0] * k) + r0->size[0] * r0->size[1] * y];
            }
        }
    }
    /* 'nondecimatedWaveletTransform3D:26' B(:, :, 6*L+1:7*L)=directionalConv(BHH, Lo_D, 3); */
    i0 = 6 * p + 1;
    if (i0 > 7 * p) {
        i0 = 0;
    } else {
        i0--;
    }
    d_directionalConv(BHH, r0);
    loop_ub = r0->size[2] - 1;
    for (y = 0; y <= loop_ub; y++) {
        for (k = 0; k < 183; k++) {
            for (b_p = 0; b_p < 263; b_p++) {
                B->data[(b_p + B->size[0] * k) + B->size[0] * B->size[1] * (i0 + y)] = r0->data[(b_p + r0->size[0] * k) + r0->size[0] * r0->size[1] * y];
            }
        }
    }
    /* 'nondecimatedWaveletTransform3D:27' B(:, :, 7*L+1:8*L)=directionalConv(BHH, Hi_D, 3); */
    i0 = 7 * p + 1;
    if (i0 > (p << 3)) {
        i0 = 0;
    } else {
        i0--;
    }
    e_directionalConv(BHH, r0);
    emxFree_creal_T(&BHH);
    loop_ub = r0->size[2] - 1;
    for (y = 0; y <= loop_ub; y++) {
        for (k = 0; k < 183; k++) {
            for (b_p = 0; b_p < 263; b_p++) {
                B->data[(b_p + B->size[0] * k) + B->size[0] * B->size[1] * (i0 + y)] = r0->data[(b_p + r0->size[0] * k) + r0->size[0] * r0->size[1] * y];
            }
        }
    }
    emxFree_creal_T(&r0);
    emlrtHeapReferenceStackLeaveFcn();
}

void nondecimatedWaveletTransform3D_api(const mxArray * const prhs[1], const mxArray *plhs[1])
{
    emxArray_creal_T *A;
    emxArray_creal_T *B;
    emlrtHeapReferenceStackEnterFcn();
    emxInit_creal_T(&A, 3, TRUE);
    emxInit_creal_T(&B, 3, TRUE);
    /* Marshall function inputs */
    emlrt_marshallIn(emlrtAliasP(prhs[0]), "A", A);
    /* Invoke the target function */
    nondecimatedWaveletTransform3D(A, B);
    /* Marshall function outputs */
    plhs[0] = emlrt_marshallOut(B);
    emxFree_creal_T(&B);
    emxFree_creal_T(&A);
    emlrtHeapReferenceStackLeaveFcn();
}

void nondecimatedWaveletTransform3D_atexit(void)
{
    emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void nondecimatedWaveletTransform3D_initialize(emlrtContext *context)
{
    emlrtEnterRtStack(&emlrtContextGlobal);
    emlrtFirstTime(context);
}

void nondecimatedWaveletTransform3D_terminate(void)
{
    emlrtLeaveRtStack(&emlrtContextGlobal);
}
/* End of code generation (nondecimatedWaveletTransform3D.c) */
