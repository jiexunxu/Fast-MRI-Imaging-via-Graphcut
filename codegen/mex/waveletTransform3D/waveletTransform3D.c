/*
 * waveletTransform3D.c
 *
 * Code generation for function 'waveletTransform3D'
 *
 * C source code generated on: Wed Sep 05 17:20:05 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "waveletTransform3D.h"

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
static emlrtDCInfo emlrtDCI = { 24, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo b_emlrtDCI = { 25, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo c_emlrtDCI = { 25, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo d_emlrtDCI = { 26, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo e_emlrtDCI = { 26, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo f_emlrtDCI = { 27, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo g_emlrtDCI = { 27, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo h_emlrtDCI = { 28, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo i_emlrtDCI = { 28, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo j_emlrtDCI = { 29, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo k_emlrtDCI = { 29, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo l_emlrtDCI = { 30, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo m_emlrtDCI = { 30, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo n_emlrtDCI = { 31, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo o_emlrtDCI = { 31, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo p_emlrtDCI = { 16, 21, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo q_emlrtDCI = { 16, 28, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo r_emlrtDCI = { 16, 35, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo s_emlrtDCI = { 16, 21, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo t_emlrtDCI = { 16, 28, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo u_emlrtDCI = { 16, 35, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo v_emlrtDCI = { 39, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo w_emlrtDCI = { 40, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo x_emlrtDCI = { 40, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo y_emlrtDCI = { 41, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo ab_emlrtDCI = { 41, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo bb_emlrtDCI = { 42, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo cb_emlrtDCI = { 42, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo db_emlrtDCI = { 43, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo eb_emlrtDCI = { 43, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo fb_emlrtDCI = { 44, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo gb_emlrtDCI = { 44, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo hb_emlrtDCI = { 45, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo ib_emlrtDCI = { 45, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo jb_emlrtDCI = { 46, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo kb_emlrtDCI = { 46, 34, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo lb_emlrtDCI = { 47, 9, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo mb_emlrtDCI = { 47, 14, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo nb_emlrtDCI = { 47, 19, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo ob_emlrtDCI = { 47, 19, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo pb_emlrtDCI = { 38, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo qb_emlrtDCI = { 38, 20, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo rb_emlrtDCI = { 38, 27, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo sb_emlrtDCI = { 38, 13, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo tb_emlrtDCI = { 38, 20, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };
static emlrtDCInfo ub_emlrtDCI = { 38, 27, "waveletTransform3D", "X:/Projects/ApproximateNullVec/deployment/utils/waveletTransform3D.m", 1 };

/* Function Declarations */
static void b_conv(const creal_T A[176], const emxArray_real_T *B, emxArray_creal_T *C);
static void b_directionalConv(const emxArray_creal_T *A, const emxArray_real_T *w, emxArray_creal_T *B);
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_creal_T *y);
static void b_emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void b_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void c_conv(const emxArray_creal_T *A, const emxArray_real_T *B, emxArray_creal_T *C);
static void c_directionalConv(const emxArray_creal_T *A, const emxArray_real_T *w, emxArray_creal_T *B);
static void c_emlrt_marshallIn(const mxArray *w1, const char_T *identifier, emxArray_real_T *y);
static void c_emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void conv(const creal_T A[256], const emxArray_real_T *B, emxArray_creal_T *C);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void directionalConv(const emxArray_creal_T *A, const emxArray_real_T *w, emxArray_creal_T *B);
static real_T e_emlrt_marshallIn(const mxArray *direction, const char_T *identifier);
static void emlrt_marshallIn(const mxArray *A, const char_T *identifier, emxArray_creal_T *y);
static const mxArray *emlrt_marshallOut(emxArray_creal_T *u);
static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize);
static void emxFree_creal_T(emxArray_creal_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
static real_T f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId);
static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_creal_T *ret);
static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static real_T i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId);
static real_T length(const emxArray_real_T *x);
static void mrdivide(const emxArray_real_T *A, real_T B, emxArray_real_T *y);
static void singleBlockRecon(const emxArray_creal_T *A, const emxArray_real_T *w1, const emxArray_real_T *w2, const emxArray_real_T *w3, emxArray_creal_T *B);
static void squeeze(const emxArray_creal_T *a, emxArray_creal_T *b);
static void waveletDecomposition3D(const emxArray_creal_T *A, const emxArray_real_T *Lo_D, const emxArray_real_T *Hi_D, emxArray_creal_T *B);
static void waveletReconstruction3D(const emxArray_creal_T *A, const emxArray_real_T *Lo_R, const emxArray_real_T *Hi_R, emxArray_creal_T *B);

/* Function Definitions */

/*
 * 
 */
static void b_conv(const creal_T A[176], const emxArray_real_T *B, emxArray_creal_T *C)
{
    int32_T nB;
    int32_T nApnB;
    int32_T k;
    int32_T jA1;
    int32_T jC;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    nB = B->size[1];
    nApnB = nB + 175;
    if (nB == 0) {
        nApnB++;
    }
    k = C->size[0] * C->size[1];
    C->size[0] = 1;
    C->size[1] = nApnB;
    emxEnsureCapacity((emxArray__common *)C, k, (int32_T)sizeof(creal_T));
    if ((B->size[1] == 0) || (nApnB == 0)) {
        k = C->size[0] * C->size[1];
        C->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)C, k, (int32_T)sizeof(creal_T));
        jA1 = C->size[1] - 1;
        for (k = 0; k <= jA1; k++) {
            C->data[C->size[0] * k].re = 0.0;
            C->data[C->size[0] * k].im = 0.0;
        }
    } else {
        for (jC = 1; jC <= nApnB; jC++) {
            if (nB < jC + 1) {
                jA1 = jC - nB;
            } else {
                jA1 = 0;
            }
            if (176 < jC) {
                jA2 = 176;
            } else {
                jA2 = jC;
            }
            s_re = 0.0;
            s_im = 0.0;
            for (k = jA1 + 1; k <= jA2; k++) {
                s_re += A[k - 1].re * B->data[jC - k];
                s_im += A[k - 1].im * B->data[jC - k];
            }
            C->data[jC - 1].re = s_re;
            C->data[jC - 1].im = s_im;
        }
    }
}

/*
 * function B=directionalConv(A, w, dir)
 */
static void b_directionalConv(const emxArray_creal_T *A, const emxArray_real_T *w, emxArray_creal_T *B)
{
    int32_T m;
    int32_T p;
    int32_T q;
    int32_T i6;
    uint32_T i;
    emxArray_creal_T *v;
    uint32_T k;
    creal_T b[176];
    int32_T i7;
    int32_T i8;
    emlrtHeapReferenceStackEnterFcn();
    /* 'waveletTransform3D:51' [m, n, p]=size(A); */
    m = A->size[0];
    p = A->size[2];
    /* 'waveletTransform3D:52' q=length(w); */
    if (w->size[1] == 0) {
        q = 0;
    } else {
        q = w->size[1];
    }
    /* 'waveletTransform3D:53' if dir==1 */
    /* 'waveletTransform3D:61' elseif dir==2 */
    /* 'waveletTransform3D:62' B=complex(zeros(m, n+q-1, p), 0); */
    i6 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = m;
    B->size[1] = (int32_T)((176.0 + (real_T)q) - 1.0);
    B->size[2] = p;
    emxEnsureCapacity((emxArray__common *)B, i6, (int32_T)sizeof(creal_T));
    q = m * (int32_T)((176.0 + (real_T)q) - 1.0) * p - 1;
    for (i6 = 0; i6 <= q; i6++) {
        B->data[i6].re = 0.0;
        B->data[i6].im = 0.0;
    }
    /* 'waveletTransform3D:63' for i=1:m */
    i = 1U;
    c_emxInit_creal_T(&v, 2, TRUE);
    while (i <= (uint32_T)m) {
        /* 'waveletTransform3D:64' for k=1:p */
        for (k = 1U; k <= (uint32_T)p; k++) {
            /* 'waveletTransform3D:65' v=conv(squeeze(A(i, :, k)), w); */
            for (q = 0; q < 176; q++) {
                b[q] = A->data[(((int32_T)i + A->size[0] * q) + A->size[0] * A->size[1] * ((int32_T)k - 1)) - 1];
            }
            b_conv(b, w, v);
            /* 'waveletTransform3D:66' B(i, :, k)=v; */
            i6 = (int32_T)i - 1;
            i7 = (int32_T)k - 1;
            q = v->size[1] - 1;
            for (i8 = 0; i8 <= q; i8++) {
                B->data[(i6 + B->size[0] * i8) + B->size[0] * B->size[1] * i7] = v->data[v->size[0] * i8];
            }
        }
        i++;
    }
    emxFree_creal_T(&v);
    emlrtHeapReferenceStackLeaveFcn();
}

static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_creal_T *y)
{
    g_emlrt_marshallIn(emlrtAlias(u), parentId, y);
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

static void b_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
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
 * 
 */
static void c_conv(const emxArray_creal_T *A, const emxArray_real_T *B, emxArray_creal_T *C)
{
    int32_T nA;
    int32_T nB;
    int32_T nApnB;
    int32_T k;
    int32_T jA1;
    int32_T jC;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    nA = A->size[0];
    nB = B->size[1];
    nApnB = nA + nB;
    if ((nA == 0) || (nB == 0)) {
    } else {
        nApnB--;
    }
    k = C->size[0];
    C->size[0] = nApnB;
    emxEnsureCapacity((emxArray__common *)C, k, (int32_T)sizeof(creal_T));
    if ((A->size[0] == 0) || (B->size[1] == 0) || (nApnB == 0)) {
        jA1 = C->size[0];
        k = C->size[0];
        C->size[0] = jA1;
        emxEnsureCapacity((emxArray__common *)C, k, (int32_T)sizeof(creal_T));
        jA1--;
        for (k = 0; k <= jA1; k++) {
            C->data[k].re = 0.0;
            C->data[k].im = 0.0;
        }
    } else {
        for (jC = 1; jC <= nApnB; jC++) {
            if (nB < jC + 1) {
                jA1 = jC - nB;
            } else {
                jA1 = 0;
            }
            if (nA < jC) {
                jA2 = nA;
            } else {
                jA2 = jC;
            }
            s_re = 0.0;
            s_im = 0.0;
            for (k = jA1 + 1; k <= jA2; k++) {
                s_re += A->data[k - 1].re * B->data[jC - k];
                s_im += A->data[k - 1].im * B->data[jC - k];
            }
            C->data[jC - 1].re = s_re;
            C->data[jC - 1].im = s_im;
        }
    }
}

/*
 * function B=directionalConv(A, w, dir)
 */
static void c_directionalConv(const emxArray_creal_T *A, const emxArray_real_T *w, emxArray_creal_T *B)
{
    int32_T m;
    int32_T n;
    int32_T p;
    int32_T q;
    int32_T i9;
    int32_T loop_ub;
    uint32_T i;
    emxArray_creal_T *v;
    emxArray_int32_T *r1;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    uint32_T j;
    int32_T iv0[3];
    int32_T b_loop_ub;
    int32_T i10;
    int32_T c_loop_ub;
    int32_T i11;
    emxArray_creal_T c_v;
    emlrtHeapReferenceStackEnterFcn();
    /* 'waveletTransform3D:51' [m, n, p]=size(A); */
    m = A->size[0];
    n = A->size[1];
    p = A->size[2];
    /* 'waveletTransform3D:52' q=length(w); */
    if (w->size[1] == 0) {
        q = 0;
    } else {
        q = w->size[1];
    }
    /* 'waveletTransform3D:53' if dir==1 */
    /* 'waveletTransform3D:69' elseif dir==3 */
    /* 'waveletTransform3D:70' B=complex(zeros(m, n, p+q-1), 0); */
    i9 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = m;
    B->size[1] = n;
    B->size[2] = (int32_T)((real_T)((uint32_T)p + (uint32_T)q) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, i9, (int32_T)sizeof(creal_T));
    loop_ub = m * n * (int32_T)((real_T)((uint32_T)p + (uint32_T)q) - 1.0) - 1;
    for (i9 = 0; i9 <= loop_ub; i9++) {
        B->data[i9].re = 0.0;
        B->data[i9].im = 0.0;
    }
    /* 'waveletTransform3D:71' for i=1:m */
    i = 1U;
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r1, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    while (i <= (uint32_T)m) {
        /* 'waveletTransform3D:72' for j=1:n */
        for (j = 1U; j <= (uint32_T)n; j++) {
            /* 'waveletTransform3D:73' v=conv(squeeze(A(i, j, :)), w); */
            p = (int32_T)j;
            q = (int32_T)i;
            i9 = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, i9, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (i9 = 0; i9 <= loop_ub; i9++) {
                b_A->data[b_A->size[0] * b_A->size[1] * i9] = A->data[((q + A->size[0] * (p - 1)) + A->size[0] * A->size[1] * i9) - 1];
            }
            squeeze(b_A, v);
            i9 = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, i9, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (i9 = 0; i9 <= loop_ub; i9++) {
                b_v->data[i9] = v->data[i9];
            }
            c_conv(b_v, w, v);
            /* 'waveletTransform3D:74' B(i, j, :)=v; */
            i9 = (int32_T)i - 1;
            p = (int32_T)j - 1;
            q = r1->size[0];
            r1->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r1, q, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (q = 0; q <= loop_ub; q++) {
                r1->data[q] = q;
            }
            iv0[0] = 1;
            iv0[1] = 1;
            iv0[2] = r1->size[0];
            loop_ub = iv0[2] - 1;
            for (q = 0; q <= loop_ub; q++) {
                b_loop_ub = iv0[1] - 1;
                for (i10 = 0; i10 <= b_loop_ub; i10++) {
                    c_loop_ub = iv0[0] - 1;
                    for (i11 = 0; i11 <= c_loop_ub; i11++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv0;
                        c_v.numDimensions = 1;
                        B->data[(i9 + B->size[0] * p) + B->size[0] * B->size[1] * r1->data[q]] = c_v.data[(i11 + c_v.size[0] * i10) + c_v.size[0] * c_v.size[1] * q];
                    }
                }
            }
        }
        i++;
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r1);
    emxFree_creal_T(&v);
    emlrtHeapReferenceStackLeaveFcn();
}

static void c_emlrt_marshallIn(const mxArray *w1, const char_T *identifier, emxArray_real_T *y)
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    d_emlrt_marshallIn(emlrtAlias(w1), &thisId, y);
    emlrtDestroyArray(&w1);
}

static void c_emxInit_creal_T(emxArray_creal_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
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
static void conv(const creal_T A[256], const emxArray_real_T *B, emxArray_creal_T *C)
{
    int32_T nB;
    int32_T nApnB;
    int32_T k;
    int32_T jA1;
    int32_T jC;
    int32_T jA2;
    real_T s_re;
    real_T s_im;
    nB = B->size[1];
    nApnB = nB + 255;
    if (nB == 0) {
        nApnB++;
    }
    k = C->size[0];
    C->size[0] = nApnB;
    emxEnsureCapacity((emxArray__common *)C, k, (int32_T)sizeof(creal_T));
    if ((B->size[1] == 0) || (nApnB == 0)) {
        jA1 = C->size[0];
        k = C->size[0];
        C->size[0] = jA1;
        emxEnsureCapacity((emxArray__common *)C, k, (int32_T)sizeof(creal_T));
        jA1--;
        for (k = 0; k <= jA1; k++) {
            C->data[k].re = 0.0;
            C->data[k].im = 0.0;
        }
    } else {
        for (jC = 1; jC <= nApnB; jC++) {
            if (nB < jC + 1) {
                jA1 = jC - nB;
            } else {
                jA1 = 0;
            }
            if (256 < jC) {
                jA2 = 256;
            } else {
                jA2 = jC;
            }
            s_re = 0.0;
            s_im = 0.0;
            for (k = jA1 + 1; k <= jA2; k++) {
                s_re += A[k - 1].re * B->data[jC - k];
                s_im += A[k - 1].im * B->data[jC - k];
            }
            C->data[jC - 1].re = s_re;
            C->data[jC - 1].im = s_im;
        }
    }
}

static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
    h_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

/*
 * function B=directionalConv(A, w, dir)
 */
static void directionalConv(const emxArray_creal_T *A, const emxArray_real_T *w, emxArray_creal_T *B)
{
    int32_T p;
    int32_T q;
    int32_T i4;
    emxArray_creal_T *v;
    int32_T j;
    uint32_T k;
    creal_T b[256];
    int32_T i5;
    emlrtHeapReferenceStackEnterFcn();
    /* 'waveletTransform3D:51' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'waveletTransform3D:52' q=length(w); */
    if (w->size[1] == 0) {
        q = 0;
    } else {
        q = w->size[1];
    }
    /* 'waveletTransform3D:53' if dir==1 */
    /* 'waveletTransform3D:54' B=complex(zeros(m+q-1, n, p), 0); */
    i4 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = (int32_T)((256.0 + (real_T)q) - 1.0);
    B->size[1] = 176;
    B->size[2] = p;
    emxEnsureCapacity((emxArray__common *)B, i4, (int32_T)sizeof(creal_T));
    q = (int32_T)((256.0 + (real_T)q) - 1.0) * 176 * p - 1;
    for (i4 = 0; i4 <= q; i4++) {
        B->data[i4].re = 0.0;
        B->data[i4].im = 0.0;
    }
    /* 'waveletTransform3D:55' for j=1:n */
    b_emxInit_creal_T(&v, 1, TRUE);
    for (j = 0; j < 176; j++) {
        /* 'waveletTransform3D:56' for k=1:p */
        for (k = 1U; k <= (uint32_T)p; k++) {
            /* 'waveletTransform3D:57' v=conv(squeeze(A(:, j, k)), w); */
            for (q = 0; q < 256; q++) {
                b[q] = A->data[(q + A->size[0] * j) + A->size[0] * A->size[1] * ((int32_T)k - 1)];
            }
            conv(b, w, v);
            /* 'waveletTransform3D:58' B(:, j, k)=v; */
            i4 = (int32_T)k - 1;
            q = v->size[0] - 1;
            for (i5 = 0; i5 <= q; i5++) {
                B->data[(i5 + B->size[0] * j) + B->size[0] * B->size[1] * i4] = v->data[i5];
            }
        }
    }
    emxFree_creal_T(&v);
    emlrtHeapReferenceStackLeaveFcn();
}

static real_T e_emlrt_marshallIn(const mxArray *direction, const char_T *identifier)
{
    real_T y;
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    y = f_emlrt_marshallIn(emlrtAlias(direction), &thisId);
    emlrtDestroyArray(&direction);
    return y;
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

static real_T f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId)
{
    real_T y;
    y = i_emlrt_marshallIn(emlrtAlias(u), parentId);
    emlrtDestroyArray(&u);
    return y;
}

static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_creal_T *ret)
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

static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
    int32_T i20;
    int32_T iv4[2];
    static const boolean_T bv2[2] = { FALSE, TRUE };
    boolean_T bv3[2];
    for (i20 = 0; i20 < 2; i20++) {
        iv4[i20] = 1 + 7 * i20;
        bv3[i20] = bv2[i20];
    }
    emlrtCheckVsBuiltInR2011a(msgId, src, "double", FALSE, 2U, iv4, bv3, ret->size);
    i20 = ret->size[0] * ret->size[1];
    ret->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)ret, i20, (int32_T)sizeof(real_T));
    emlrtImportArrayR2008b(src, ret->data, 8);
    emlrtDestroyArray(&src);
}

static real_T i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId)
{
    real_T ret;
    emlrtCheckBuiltInR2011a(msgId, src, "double", FALSE, 0U, 0);
    ret = *(real_T *)mxGetData(src);
    emlrtDestroyArray(&src);
    return ret;
}

/*
 * 
 */
static real_T length(const emxArray_real_T *x)
{
    real_T n;
    if (x->size[1] == 0) {
        n = 0.0;
    } else {
        n = (real_T)x->size[1];
    }
    return n;
}

/*
 * 
 */
static void mrdivide(const emxArray_real_T *A, real_T B, emxArray_real_T *y)
{
    int32_T i0;
    int32_T loop_ub;
    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    loop_ub = A->size[0] * A->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        y->data[i0] = A->data[i0] / B;
    }
}

/*
 * function B=singleBlockRecon(A, w1, w2, w3)
 */
static void singleBlockRecon(const emxArray_creal_T *A, const emxArray_real_T *w1, const emxArray_real_T *w2, const emxArray_real_T *w3, emxArray_creal_T *B)
{
    int32_T p;
    int32_T q;
    int32_T k;
    int32_T loop_ub;
    emxArray_creal_T *v;
    emxArray_int32_T *r7;
    emxArray_creal_T *b_v;
    emxArray_creal_T *b_A;
    int32_T i;
    int32_T j;
    int32_T iv1[3];
    int32_T b_loop_ub;
    int32_T i18;
    int32_T c_loop_ub;
    int32_T i19;
    emxArray_creal_T c_v;
    emxArray_creal_T *d_v;
    real_T d_loop_ub;
    uint32_T b_k;
    creal_T b[176];
    uint32_T e_loop_ub;
    uint32_T b_j;
    creal_T b_b[256];
    emlrtHeapReferenceStackEnterFcn();
    /*  Conv in z, y, x direction with w1, w2, w3 */
    /* 'waveletTransform3D:82' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'waveletTransform3D:83' q=length(w1); */
    if (w1->size[1] == 0) {
        q = 0;
    } else {
        q = w1->size[1];
    }
    /* 'waveletTransform3D:84' B=complex(zeros(m+q-1, n+q-1, p+q-1), 0); */
    k = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = (int32_T)((256.0 + (real_T)q) - 1.0);
    B->size[1] = (int32_T)((176.0 + (real_T)q) - 1.0);
    B->size[2] = (int32_T)((real_T)((uint32_T)p + (uint32_T)q) - 1.0);
    emxEnsureCapacity((emxArray__common *)B, k, (int32_T)sizeof(creal_T));
    loop_ub = (int32_T)((256.0 + (real_T)q) - 1.0) * (int32_T)((176.0 + (real_T)q) - 1.0) * (int32_T)((real_T)((uint32_T)p + (uint32_T)q) - 1.0) - 1;
    for (k = 0; k <= loop_ub; k++) {
        B->data[k].re = 0.0;
        B->data[k].im = 0.0;
    }
    /* 'waveletTransform3D:85' for i=1:m */
    b_emxInit_creal_T(&v, 1, TRUE);
    emxInit_int32_T(&r7, 1, TRUE);
    b_emxInit_creal_T(&b_v, 1, TRUE);
    emxInit_creal_T(&b_A, 3, TRUE);
    for (i = 0; i < 256; i++) {
        /* 'waveletTransform3D:86' for j=1:n */
        for (j = 0; j < 176; j++) {
            /* 'waveletTransform3D:87' v=conv(squeeze(A(i, j, :)), w1); */
            k = b_A->size[0] * b_A->size[1] * b_A->size[2];
            b_A->size[0] = 1;
            b_A->size[1] = 1;
            b_A->size[2] = A->size[2];
            emxEnsureCapacity((emxArray__common *)b_A, k, (int32_T)sizeof(creal_T));
            loop_ub = A->size[2] - 1;
            for (k = 0; k <= loop_ub; k++) {
                b_A->data[b_A->size[0] * b_A->size[1] * k] = A->data[(i + A->size[0] * j) + A->size[0] * A->size[1] * k];
            }
            squeeze(b_A, v);
            k = b_v->size[0];
            b_v->size[0] = v->size[0];
            emxEnsureCapacity((emxArray__common *)b_v, k, (int32_T)sizeof(creal_T));
            loop_ub = v->size[0] - 1;
            for (k = 0; k <= loop_ub; k++) {
                b_v->data[k] = v->data[k];
            }
            c_conv(b_v, w1, v);
            /* 'waveletTransform3D:88' B(i, j, :)=v; */
            k = r7->size[0];
            r7->size[0] = B->size[2];
            emxEnsureCapacity((emxArray__common *)r7, k, (int32_T)sizeof(int32_T));
            loop_ub = B->size[2] - 1;
            for (k = 0; k <= loop_ub; k++) {
                r7->data[k] = k;
            }
            iv1[0] = 1;
            iv1[1] = 1;
            iv1[2] = r7->size[0];
            loop_ub = iv1[2] - 1;
            for (k = 0; k <= loop_ub; k++) {
                b_loop_ub = iv1[1] - 1;
                for (i18 = 0; i18 <= b_loop_ub; i18++) {
                    c_loop_ub = iv1[0] - 1;
                    for (i19 = 0; i19 <= c_loop_ub; i19++) {
                        c_v = *v;
                        c_v.size = (int32_T *)&iv1;
                        c_v.numDimensions = 1;
                        B->data[(i + B->size[0] * j) + B->size[0] * B->size[1] * r7->data[k]] = c_v.data[(i19 + c_v.size[0] * i18) + c_v.size[0] * c_v.size[1] * k];
                    }
                }
            }
        }
    }
    emxFree_creal_T(&b_A);
    emxFree_creal_T(&b_v);
    emxFree_int32_T(&r7);
    /* 'waveletTransform3D:92' for i=1:m */
    c_emxInit_creal_T(&d_v, 2, TRUE);
    for (i = 0; i < 256; i++) {
        /* 'waveletTransform3D:93' for k=1:p+q-1 */
        d_loop_ub = (real_T)((uint32_T)p + (uint32_T)q) - 1.0;
        for (b_k = 1U; (real_T)b_k <= d_loop_ub; b_k++) {
            /* 'waveletTransform3D:94' v=conv(squeeze(B(i, 1:n, k)), w2); */
            for (k = 0; k < 176; k++) {
                b[k] = B->data[(i + B->size[0] * k) + B->size[0] * B->size[1] * ((int32_T)b_k - 1)];
            }
            b_conv(b, w2, d_v);
            /* 'waveletTransform3D:95' B(i, :, k)=v; */
            k = (int32_T)b_k - 1;
            loop_ub = d_v->size[1] - 1;
            for (i18 = 0; i18 <= loop_ub; i18++) {
                B->data[(i + B->size[0] * i18) + B->size[0] * B->size[1] * k] = d_v->data[d_v->size[0] * i18];
            }
        }
    }
    emxFree_creal_T(&d_v);
    /* 'waveletTransform3D:99' for j=1:n+q-1 */
    e_loop_ub = (uint32_T)q + 175U;
    for (b_j = 1U; b_j <= e_loop_ub; b_j++) {
        /* 'waveletTransform3D:100' for k=1:p+q-1 */
        d_loop_ub = (real_T)((uint32_T)p + (uint32_T)q) - 1.0;
        for (b_k = 1U; (real_T)b_k <= d_loop_ub; b_k++) {
            /* 'waveletTransform3D:101' v=conv(squeeze(B(1:m, j, k)), w3); */
            for (k = 0; k < 256; k++) {
                b_b[k] = B->data[(k + B->size[0] * ((int32_T)b_j - 1)) + B->size[0] * B->size[1] * ((int32_T)b_k - 1)];
            }
            conv(b_b, w3, v);
            /* 'waveletTransform3D:102' B(:, j, k)=v; */
            k = (int32_T)b_j - 1;
            i18 = (int32_T)b_k - 1;
            loop_ub = v->size[0] - 1;
            for (i19 = 0; i19 <= loop_ub; i19++) {
                B->data[(i19 + B->size[0] * k) + B->size[0] * B->size[1] * i18] = v->data[i19];
            }
        }
    }
    emxFree_creal_T(&v);
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
 * function B=waveletDecomposition3D(A, Lo_D, Hi_D)
 */
static void waveletDecomposition3D(const emxArray_creal_T *A, const emxArray_real_T *Lo_D, const emxArray_real_T *Hi_D, emxArray_creal_T *B)
{
    int32_T p;
    real_T q;
    real_T y;
    int32_T i1;
    int32_T loop_ub;
    emxArray_creal_T *BL;
    emxArray_creal_T *BH;
    emxArray_creal_T *BLL;
    emxArray_creal_T *BLH;
    emxArray_creal_T *BHL;
    emxArray_creal_T *BHH;
    emxArray_creal_T *r0;
    int32_T b_loop_ub;
    int32_T c_loop_ub;
    int32_T i2;
    real_T d0;
    int32_T i3;
    emlrtHeapReferenceStackEnterFcn();
    /*  LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH */
    /* 'waveletTransform3D:14' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'waveletTransform3D:15' q=length(Lo_D); */
    q = length(Lo_D);
    /* 'waveletTransform3D:16' B=complex(zeros(m+q-1, n+q-1, (p+q-1)*8), 0); */
    y = (((real_T)p + q) - 1.0) * 8.0;
    i1 = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = (int32_T)emlrtIntegerCheckR2011a((256.0 + q) - 1.0, &p_emlrtDCI, &emlrtContextGlobal);
    B->size[1] = (int32_T)emlrtIntegerCheckR2011a((176.0 + q) - 1.0, &q_emlrtDCI, &emlrtContextGlobal);
    B->size[2] = (int32_T)emlrtIntegerCheckR2011a(y, &r_emlrtDCI, &emlrtContextGlobal);
    emxEnsureCapacity((emxArray__common *)B, i1, (int32_T)sizeof(creal_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((256.0 + q) - 1.0, &s_emlrtDCI, &emlrtContextGlobal) * (int32_T)emlrtIntegerCheckR2011a((176.0 + q) - 1.0, &t_emlrtDCI, &emlrtContextGlobal) * (int32_T)emlrtIntegerCheckR2011a(y, &u_emlrtDCI, &emlrtContextGlobal) - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
        B->data[i1].re = 0.0;
        B->data[i1].im = 0.0;
    }
    emxInit_creal_T(&BL, 3, TRUE);
    emxInit_creal_T(&BH, 3, TRUE);
    emxInit_creal_T(&BLL, 3, TRUE);
    emxInit_creal_T(&BLH, 3, TRUE);
    emxInit_creal_T(&BHL, 3, TRUE);
    emxInit_creal_T(&BHH, 3, TRUE);
    /* 'waveletTransform3D:17' BL=directionalConv(A, Lo_D, 1); */
    directionalConv(A, Lo_D, BL);
    /* 'waveletTransform3D:18' BH=directionalConv(A, Hi_D, 1); */
    directionalConv(A, Hi_D, BH);
    /* 'waveletTransform3D:19' BLL=directionalConv(BL, Lo_D, 2); */
    b_directionalConv(BL, Lo_D, BLL);
    /* 'waveletTransform3D:20' BLH=directionalConv(BL, Hi_D, 2); */
    b_directionalConv(BL, Hi_D, BLH);
    /* 'waveletTransform3D:21' BHL=directionalConv(BH, Lo_D, 2); */
    b_directionalConv(BH, Lo_D, BHL);
    /* 'waveletTransform3D:22' BHH=directionalConv(BH, Hi_D, 2); */
    b_directionalConv(BH, Hi_D, BHH);
    /* 'waveletTransform3D:23' L=p+q-1; */
    q = ((real_T)p + q) - 1.0;
    /* 'waveletTransform3D:24' B(:, :, 1:L)=directionalConv(BLL, Lo_D, 3); */
    emxFree_creal_T(&BH);
    emxFree_creal_T(&BL);
    if (1.0 > q) {
    } else {
        emlrtIntegerCheckR2011a(q, &emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&r0, 3, TRUE);
    c_directionalConv(BLL, Lo_D, r0);
    loop_ub = r0->size[2] - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
        b_loop_ub = r0->size[1] - 1;
        for (p = 0; p <= b_loop_ub; p++) {
            c_loop_ub = r0->size[0] - 1;
            for (i2 = 0; i2 <= c_loop_ub; i2++) {
                B->data[(i2 + B->size[0] * p) + B->size[0] * B->size[1] * i1] = r0->data[(i2 + r0->size[0] * p) + r0->size[0] * r0->size[1] * i1];
            }
        }
    }
    /* 'waveletTransform3D:25' B(:, :, L+1:2*L)=directionalConv(BLL, Hi_D, 3); */
    y = q + 1.0;
    d0 = 2.0 * q;
    if (y > d0) {
        i1 = 0;
    } else {
        i1 = (int32_T)emlrtIntegerCheckR2011a(y, &b_emlrtDCI, &emlrtContextGlobal) - 1;
        emlrtIntegerCheckR2011a(d0, &c_emlrtDCI, &emlrtContextGlobal);
    }
    c_directionalConv(BLL, Hi_D, r0);
    emxFree_creal_T(&BLL);
    loop_ub = r0->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_loop_ub = r0->size[1] - 1;
        for (i2 = 0; i2 <= b_loop_ub; i2++) {
            c_loop_ub = r0->size[0] - 1;
            for (i3 = 0; i3 <= c_loop_ub; i3++) {
                B->data[(i3 + B->size[0] * i2) + B->size[0] * B->size[1] * (i1 + p)] = r0->data[(i3 + r0->size[0] * i2) + r0->size[0] * r0->size[1] * p];
            }
        }
    }
    /* 'waveletTransform3D:26' B(:, :, 2*L+1:3*L)=directionalConv(BLH, Lo_D, 3); */
    y = 2.0 * q + 1.0;
    d0 = 3.0 * q;
    if (y > d0) {
        i1 = 0;
    } else {
        i1 = (int32_T)emlrtIntegerCheckR2011a(y, &d_emlrtDCI, &emlrtContextGlobal) - 1;
        emlrtIntegerCheckR2011a(d0, &e_emlrtDCI, &emlrtContextGlobal);
    }
    c_directionalConv(BLH, Lo_D, r0);
    loop_ub = r0->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_loop_ub = r0->size[1] - 1;
        for (i2 = 0; i2 <= b_loop_ub; i2++) {
            c_loop_ub = r0->size[0] - 1;
            for (i3 = 0; i3 <= c_loop_ub; i3++) {
                B->data[(i3 + B->size[0] * i2) + B->size[0] * B->size[1] * (i1 + p)] = r0->data[(i3 + r0->size[0] * i2) + r0->size[0] * r0->size[1] * p];
            }
        }
    }
    /* 'waveletTransform3D:27' B(:, :, 3*L+1:4*L)=directionalConv(BLH, Hi_D, 3); */
    y = 3.0 * q + 1.0;
    d0 = 4.0 * q;
    if (y > d0) {
        i1 = 0;
    } else {
        i1 = (int32_T)emlrtIntegerCheckR2011a(y, &f_emlrtDCI, &emlrtContextGlobal) - 1;
        emlrtIntegerCheckR2011a(d0, &g_emlrtDCI, &emlrtContextGlobal);
    }
    c_directionalConv(BLH, Hi_D, r0);
    emxFree_creal_T(&BLH);
    loop_ub = r0->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_loop_ub = r0->size[1] - 1;
        for (i2 = 0; i2 <= b_loop_ub; i2++) {
            c_loop_ub = r0->size[0] - 1;
            for (i3 = 0; i3 <= c_loop_ub; i3++) {
                B->data[(i3 + B->size[0] * i2) + B->size[0] * B->size[1] * (i1 + p)] = r0->data[(i3 + r0->size[0] * i2) + r0->size[0] * r0->size[1] * p];
            }
        }
    }
    /* 'waveletTransform3D:28' B(:, :, 4*L+1:5*L)=directionalConv(BHL, Lo_D, 3); */
    y = 4.0 * q + 1.0;
    d0 = 5.0 * q;
    if (y > d0) {
        i1 = 0;
    } else {
        i1 = (int32_T)emlrtIntegerCheckR2011a(y, &h_emlrtDCI, &emlrtContextGlobal) - 1;
        emlrtIntegerCheckR2011a(d0, &i_emlrtDCI, &emlrtContextGlobal);
    }
    c_directionalConv(BHL, Lo_D, r0);
    loop_ub = r0->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_loop_ub = r0->size[1] - 1;
        for (i2 = 0; i2 <= b_loop_ub; i2++) {
            c_loop_ub = r0->size[0] - 1;
            for (i3 = 0; i3 <= c_loop_ub; i3++) {
                B->data[(i3 + B->size[0] * i2) + B->size[0] * B->size[1] * (i1 + p)] = r0->data[(i3 + r0->size[0] * i2) + r0->size[0] * r0->size[1] * p];
            }
        }
    }
    /* 'waveletTransform3D:29' B(:, :, 5*L+1:6*L)=directionalConv(BHL, Hi_D, 3); */
    y = 5.0 * q + 1.0;
    d0 = 6.0 * q;
    if (y > d0) {
        i1 = 0;
    } else {
        i1 = (int32_T)emlrtIntegerCheckR2011a(y, &j_emlrtDCI, &emlrtContextGlobal) - 1;
        emlrtIntegerCheckR2011a(d0, &k_emlrtDCI, &emlrtContextGlobal);
    }
    c_directionalConv(BHL, Hi_D, r0);
    emxFree_creal_T(&BHL);
    loop_ub = r0->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_loop_ub = r0->size[1] - 1;
        for (i2 = 0; i2 <= b_loop_ub; i2++) {
            c_loop_ub = r0->size[0] - 1;
            for (i3 = 0; i3 <= c_loop_ub; i3++) {
                B->data[(i3 + B->size[0] * i2) + B->size[0] * B->size[1] * (i1 + p)] = r0->data[(i3 + r0->size[0] * i2) + r0->size[0] * r0->size[1] * p];
            }
        }
    }
    /* 'waveletTransform3D:30' B(:, :, 6*L+1:7*L)=directionalConv(BHH, Lo_D, 3); */
    y = 6.0 * q + 1.0;
    d0 = 7.0 * q;
    if (y > d0) {
        i1 = 0;
    } else {
        i1 = (int32_T)emlrtIntegerCheckR2011a(y, &l_emlrtDCI, &emlrtContextGlobal) - 1;
        emlrtIntegerCheckR2011a(d0, &m_emlrtDCI, &emlrtContextGlobal);
    }
    c_directionalConv(BHH, Lo_D, r0);
    loop_ub = r0->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_loop_ub = r0->size[1] - 1;
        for (i2 = 0; i2 <= b_loop_ub; i2++) {
            c_loop_ub = r0->size[0] - 1;
            for (i3 = 0; i3 <= c_loop_ub; i3++) {
                B->data[(i3 + B->size[0] * i2) + B->size[0] * B->size[1] * (i1 + p)] = r0->data[(i3 + r0->size[0] * i2) + r0->size[0] * r0->size[1] * p];
            }
        }
    }
    /* 'waveletTransform3D:31' B(:, :, 7*L+1:8*L)=directionalConv(BHH, Hi_D, 3); */
    y = 7.0 * q + 1.0;
    d0 = 8.0 * q;
    if (y > d0) {
        i1 = 0;
    } else {
        i1 = (int32_T)emlrtIntegerCheckR2011a(y, &n_emlrtDCI, &emlrtContextGlobal) - 1;
        emlrtIntegerCheckR2011a(d0, &o_emlrtDCI, &emlrtContextGlobal);
    }
    c_directionalConv(BHH, Hi_D, r0);
    emxFree_creal_T(&BHH);
    loop_ub = r0->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_loop_ub = r0->size[1] - 1;
        for (i2 = 0; i2 <= b_loop_ub; i2++) {
            c_loop_ub = r0->size[0] - 1;
            for (i3 = 0; i3 <= c_loop_ub; i3++) {
                B->data[(i3 + B->size[0] * i2) + B->size[0] * B->size[1] * (i1 + p)] = r0->data[(i3 + r0->size[0] * i2) + r0->size[0] * r0->size[1] * p];
            }
        }
    }
    emxFree_creal_T(&r0);
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * function B=waveletReconstruction3D(A, Lo_R, Hi_R)
 */
static void waveletReconstruction3D(const emxArray_creal_T *A, const emxArray_real_T *Lo_R, const emxArray_real_T *Hi_R, emxArray_creal_T *B)
{
    emxArray_real_T *b_B;
    int32_T p;
    real_T b_p;
    real_T q;
    real_T L;
    int32_T loop_ub;
    emxArray_creal_T *b_A;
    int32_T c_B;
    int32_T i12;
    real_T d1;
    real_T d2;
    emxArray_creal_T *c_A;
    int32_T i13;
    emxArray_creal_T *r2;
    emxArray_creal_T *d_A;
    emxArray_creal_T *r3;
    emxArray_creal_T *e_A;
    emxArray_creal_T *r4;
    emxArray_creal_T *f_A;
    emxArray_creal_T *r5;
    emxArray_creal_T *g_A;
    emxArray_creal_T *r6;
    emxArray_creal_T *h_A;
    emxArray_creal_T *i_A;
    int32_T i14;
    int32_T i15;
    emxArray_creal_T *d_B;
    int32_T i16;
    int32_T b_loop_ub;
    int32_T c_loop_ub;
    int32_T i17;
    emlrtHeapReferenceStackEnterFcn();
    b_emxInit_real_T(&b_B, 3, TRUE);
    /* 'waveletTransform3D:35' [m, n, p]=size(A); */
    p = A->size[2];
    /* 'waveletTransform3D:35' p=p/8; */
    b_p = (real_T)p / 8.0;
    /* 'waveletTransform3D:36' q=length(Lo_R); */
    q = length(Lo_R);
    /* 'waveletTransform3D:37' L=p+q-1; */
    L = (b_p + q) - 1.0;
    /* 'waveletTransform3D:38' B=zeros(m+q-1, n+q-1, p+q-1); */
    p = b_B->size[0] * b_B->size[1] * b_B->size[2];
    b_B->size[0] = (int32_T)emlrtIntegerCheckR2011a((256.0 + q) - 1.0, &pb_emlrtDCI, &emlrtContextGlobal);
    b_B->size[1] = (int32_T)emlrtIntegerCheckR2011a((176.0 + q) - 1.0, &qb_emlrtDCI, &emlrtContextGlobal);
    b_B->size[2] = (int32_T)emlrtIntegerCheckR2011a((b_p + q) - 1.0, &rb_emlrtDCI, &emlrtContextGlobal);
    emxEnsureCapacity((emxArray__common *)b_B, p, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((256.0 + q) - 1.0, &sb_emlrtDCI, &emlrtContextGlobal) * (int32_T)emlrtIntegerCheckR2011a((176.0 + q) - 1.0, &tb_emlrtDCI, &emlrtContextGlobal) * (int32_T)emlrtIntegerCheckR2011a((b_p + q) - 1.0, &ub_emlrtDCI, &emlrtContextGlobal) - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_B->data[p] = 0.0;
    }
    /* 'waveletTransform3D:39' B=B+singleBlockRecon(A(:, :, 1:L), Lo_R, Lo_R, Lo_R); */
    if (1.0 > L) {
        p = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(L, &v_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&b_A, 3, TRUE);
    c_B = b_A->size[0] * b_A->size[1] * b_A->size[2];
    b_A->size[0] = 256;
    b_A->size[1] = 176;
    b_A->size[2] = p;
    emxEnsureCapacity((emxArray__common *)b_A, c_B, (int32_T)sizeof(creal_T));
    loop_ub = p - 1;
    for (p = 0; p <= loop_ub; p++) {
        for (c_B = 0; c_B < 176; c_B++) {
            for (i12 = 0; i12 < 256; i12++) {
                b_A->data[(i12 + b_A->size[0] * c_B) + b_A->size[0] * b_A->size[1] * p] = A->data[(i12 + A->size[0] * c_B) + A->size[0] * A->size[1] * p];
            }
        }
    }
    singleBlockRecon(b_A, Lo_R, Lo_R, Lo_R, B);
    /* 'waveletTransform3D:40' B=B+singleBlockRecon(A(:, :, L+1:2*L), Hi_R, Lo_R, Lo_R); */
    d1 = L + 1.0;
    d2 = 2.0 * L;
    emxFree_creal_T(&b_A);
    if (d1 > d2) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d1, &w_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d2, &x_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&c_A, 3, TRUE);
    i12 = c_A->size[0] * c_A->size[1] * c_A->size[2];
    c_A->size[0] = 256;
    c_A->size[1] = 176;
    c_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)c_A, i12, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i12 = 0; i12 < 176; i12++) {
            for (i13 = 0; i13 < 256; i13++) {
                c_A->data[(i13 + c_A->size[0] * i12) + c_A->size[0] * c_A->size[1] * c_B] = A->data[(i13 + A->size[0] * i12) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r2, 3, TRUE);
    singleBlockRecon(c_A, Hi_R, Lo_R, Lo_R, r2);
    /* 'waveletTransform3D:41' B=B+singleBlockRecon(A(:, :, 2*L+1:3*L), Lo_R, Hi_R, Lo_R); */
    d1 = 2.0 * L + 1.0;
    d2 = 3.0 * L;
    emxFree_creal_T(&c_A);
    if (d1 > d2) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d1, &y_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d2, &ab_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&d_A, 3, TRUE);
    i12 = d_A->size[0] * d_A->size[1] * d_A->size[2];
    d_A->size[0] = 256;
    d_A->size[1] = 176;
    d_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)d_A, i12, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i12 = 0; i12 < 176; i12++) {
            for (i13 = 0; i13 < 256; i13++) {
                d_A->data[(i13 + d_A->size[0] * i12) + d_A->size[0] * d_A->size[1] * c_B] = A->data[(i13 + A->size[0] * i12) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r3, 3, TRUE);
    singleBlockRecon(d_A, Lo_R, Hi_R, Lo_R, r3);
    /* 'waveletTransform3D:42' B=B+singleBlockRecon(A(:, :, 3*L+1:4*L), Hi_R, Hi_R, Lo_R); */
    d1 = 3.0 * L + 1.0;
    d2 = 4.0 * L;
    emxFree_creal_T(&d_A);
    if (d1 > d2) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d1, &bb_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d2, &cb_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&e_A, 3, TRUE);
    i12 = e_A->size[0] * e_A->size[1] * e_A->size[2];
    e_A->size[0] = 256;
    e_A->size[1] = 176;
    e_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)e_A, i12, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i12 = 0; i12 < 176; i12++) {
            for (i13 = 0; i13 < 256; i13++) {
                e_A->data[(i13 + e_A->size[0] * i12) + e_A->size[0] * e_A->size[1] * c_B] = A->data[(i13 + A->size[0] * i12) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r4, 3, TRUE);
    singleBlockRecon(e_A, Hi_R, Hi_R, Lo_R, r4);
    /* 'waveletTransform3D:43' B=B+singleBlockRecon(A(:, :, 4*L+1:5*L), Lo_R, Lo_R, Hi_R); */
    d1 = 4.0 * L + 1.0;
    d2 = 5.0 * L;
    emxFree_creal_T(&e_A);
    if (d1 > d2) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d1, &db_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d2, &eb_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&f_A, 3, TRUE);
    i12 = f_A->size[0] * f_A->size[1] * f_A->size[2];
    f_A->size[0] = 256;
    f_A->size[1] = 176;
    f_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)f_A, i12, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i12 = 0; i12 < 176; i12++) {
            for (i13 = 0; i13 < 256; i13++) {
                f_A->data[(i13 + f_A->size[0] * i12) + f_A->size[0] * f_A->size[1] * c_B] = A->data[(i13 + A->size[0] * i12) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r5, 3, TRUE);
    singleBlockRecon(f_A, Lo_R, Lo_R, Hi_R, r5);
    /* 'waveletTransform3D:44' B=B+singleBlockRecon(A(:, :, 5*L+1:6*L), Hi_R, Lo_R, Hi_R); */
    d1 = 5.0 * L + 1.0;
    d2 = 6.0 * L;
    emxFree_creal_T(&f_A);
    if (d1 > d2) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d1, &fb_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d2, &gb_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&g_A, 3, TRUE);
    i12 = g_A->size[0] * g_A->size[1] * g_A->size[2];
    g_A->size[0] = 256;
    g_A->size[1] = 176;
    g_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)g_A, i12, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i12 = 0; i12 < 176; i12++) {
            for (i13 = 0; i13 < 256; i13++) {
                g_A->data[(i13 + g_A->size[0] * i12) + g_A->size[0] * g_A->size[1] * c_B] = A->data[(i13 + A->size[0] * i12) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    emxInit_creal_T(&r6, 3, TRUE);
    singleBlockRecon(g_A, Hi_R, Lo_R, Hi_R, r6);
    p = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = b_B->size[0];
    B->size[1] = b_B->size[1];
    B->size[2] = b_B->size[2];
    emxEnsureCapacity((emxArray__common *)B, p, (int32_T)sizeof(creal_T));
    emxFree_creal_T(&g_A);
    loop_ub = b_B->size[0] * b_B->size[1] * b_B->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        B->data[p].re = (((((b_B->data[p] + B->data[p].re) + r2->data[p].re) + r3->data[p].re) + r4->data[p].re) + r5->data[p].re) + r6->data[p].re;
        B->data[p].im = ((((B->data[p].im + r2->data[p].im) + r3->data[p].im) + r4->data[p].im) + r5->data[p].im) + r6->data[p].im;
    }
    emxFree_creal_T(&r6);
    emxFree_creal_T(&r5);
    emxFree_creal_T(&r4);
    emxFree_creal_T(&r3);
    emxFree_real_T(&b_B);
    /* 'waveletTransform3D:45' B=B+singleBlockRecon(A(:, :, 6*L+1:7*L), Lo_R, Hi_R, Hi_R); */
    d1 = 6.0 * L + 1.0;
    d2 = 7.0 * L;
    if (d1 > d2) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d1, &hb_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d2, &ib_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&h_A, 3, TRUE);
    i12 = h_A->size[0] * h_A->size[1] * h_A->size[2];
    h_A->size[0] = 256;
    h_A->size[1] = 176;
    h_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)h_A, i12, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i12 = 0; i12 < 176; i12++) {
            for (i13 = 0; i13 < 256; i13++) {
                h_A->data[(i13 + h_A->size[0] * i12) + h_A->size[0] * h_A->size[1] * c_B] = A->data[(i13 + A->size[0] * i12) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    singleBlockRecon(h_A, Lo_R, Hi_R, Hi_R, r2);
    p = B->size[0] * B->size[1] * B->size[2];
    emxEnsureCapacity((emxArray__common *)B, p, (int32_T)sizeof(creal_T));
    p = B->size[0];
    loop_ub = B->size[1];
    c_B = B->size[2];
    emxFree_creal_T(&h_A);
    loop_ub = p * loop_ub * c_B - 1;
    for (p = 0; p <= loop_ub; p++) {
        B->data[p].re += r2->data[p].re;
        B->data[p].im += r2->data[p].im;
    }
    /* 'waveletTransform3D:46' B=B+singleBlockRecon(A(:, :, 7*L+1:8*L), Hi_R, Hi_R, Hi_R); */
    d1 = 7.0 * L + 1.0;
    d2 = 8.0 * L;
    if (d1 > d2) {
        p = 0;
        c_B = 0;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(d1, &jb_emlrtDCI, &emlrtContextGlobal) - 1;
        c_B = (int32_T)emlrtIntegerCheckR2011a(d2, &kb_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&i_A, 3, TRUE);
    i12 = i_A->size[0] * i_A->size[1] * i_A->size[2];
    i_A->size[0] = 256;
    i_A->size[1] = 176;
    i_A->size[2] = c_B - p;
    emxEnsureCapacity((emxArray__common *)i_A, i12, (int32_T)sizeof(creal_T));
    loop_ub = (c_B - p) - 1;
    for (c_B = 0; c_B <= loop_ub; c_B++) {
        for (i12 = 0; i12 < 176; i12++) {
            for (i13 = 0; i13 < 256; i13++) {
                i_A->data[(i13 + i_A->size[0] * i12) + i_A->size[0] * i_A->size[1] * c_B] = A->data[(i13 + A->size[0] * i12) + A->size[0] * A->size[1] * (p + c_B)];
            }
        }
    }
    singleBlockRecon(i_A, Hi_R, Hi_R, Hi_R, r2);
    p = B->size[0] * B->size[1] * B->size[2];
    emxEnsureCapacity((emxArray__common *)B, p, (int32_T)sizeof(creal_T));
    p = B->size[0];
    loop_ub = B->size[1];
    c_B = B->size[2];
    emxFree_creal_T(&i_A);
    loop_ub = p * loop_ub * c_B - 1;
    for (p = 0; p <= loop_ub; p++) {
        B->data[p].re += r2->data[p].re;
        B->data[p].im += r2->data[p].im;
    }
    emxFree_creal_T(&r2);
    /* 'waveletTransform3D:47' B=B(q:m, q:n, q:p); */
    if (q > 256.0) {
        p = 1;
        c_B = 1;
    } else {
        p = (int32_T)emlrtIntegerCheckR2011a(q, &lb_emlrtDCI, &emlrtContextGlobal);
        c_B = 257;
    }
    if (q > 176.0) {
        i12 = 1;
        i13 = 1;
    } else {
        i12 = (int32_T)emlrtIntegerCheckR2011a(q, &mb_emlrtDCI, &emlrtContextGlobal);
        i13 = 177;
    }
    if (q > b_p) {
        i14 = 0;
        i15 = 0;
    } else {
        i14 = (int32_T)emlrtIntegerCheckR2011a(q, &nb_emlrtDCI, &emlrtContextGlobal) - 1;
        i15 = (int32_T)emlrtIntegerCheckR2011a(b_p, &ob_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_creal_T(&d_B, 3, TRUE);
    i16 = d_B->size[0] * d_B->size[1] * d_B->size[2];
    d_B->size[0] = c_B - p;
    d_B->size[1] = i13 - i12;
    d_B->size[2] = i15 - i14;
    emxEnsureCapacity((emxArray__common *)d_B, i16, (int32_T)sizeof(creal_T));
    loop_ub = (i15 - i14) - 1;
    for (i15 = 0; i15 <= loop_ub; i15++) {
        b_loop_ub = (i13 - i12) - 1;
        for (i16 = 0; i16 <= b_loop_ub; i16++) {
            c_loop_ub = (c_B - p) - 1;
            for (i17 = 0; i17 <= c_loop_ub; i17++) {
                d_B->data[(i17 + d_B->size[0] * i16) + d_B->size[0] * d_B->size[1] * i15] = B->data[(((p + i17) + B->size[0] * ((i12 + i16) - 1)) + B->size[0] * B->size[1] * (i14 + i15)) - 1];
            }
        }
    }
    p = B->size[0] * B->size[1] * B->size[2];
    B->size[0] = d_B->size[0];
    B->size[1] = d_B->size[1];
    B->size[2] = d_B->size[2];
    emxEnsureCapacity((emxArray__common *)B, p, (int32_T)sizeof(creal_T));
    loop_ub = d_B->size[2] - 1;
    for (p = 0; p <= loop_ub; p++) {
        b_loop_ub = d_B->size[1] - 1;
        for (c_B = 0; c_B <= b_loop_ub; c_B++) {
            c_loop_ub = d_B->size[0] - 1;
            for (i12 = 0; i12 <= c_loop_ub; i12++) {
                B->data[(i12 + B->size[0] * c_B) + B->size[0] * B->size[1] * p] = d_B->data[(i12 + d_B->size[0] * c_B) + d_B->size[0] * d_B->size[1] * p];
            }
        }
    }
    emxFree_creal_T(&d_B);
    emlrtHeapReferenceStackLeaveFcn();
}

/*
 * function B=waveletTransform3D(A, w1, w2, direction)
 */
void waveletTransform3D(const emxArray_creal_T *A, const emxArray_real_T *w1, const emxArray_real_T *w2, real_T direction, emxArray_creal_T *B)
{
    real_T alpha;
    int32_T vlen;
    int32_T k;
    emxArray_real_T *b_w1;
    emxArray_real_T *b_w2;
    emlrtHeapReferenceStackEnterFcn();
    /* 'waveletTransform3D:2' alpha=sum(w1); */
    if (w1->size[1] == 0) {
        alpha = 0.0;
    } else {
        vlen = w1->size[1];
        alpha = w1->data[0];
        for (k = 2; k <= vlen; k++) {
            alpha += w1->data[k - 1];
        }
    }
    emxInit_real_T(&b_w1, 2, TRUE);
    emxInit_real_T(&b_w2, 2, TRUE);
    /* 'waveletTransform3D:2' w1=w1/alpha; */
    mrdivide(w1, alpha, b_w1);
    /* 'waveletTransform3D:2' w2=w2/alpha; */
    mrdivide(w2, alpha, b_w2);
    /* 'waveletTransform3D:3' if direction==+1 */
    if (direction == 1.0) {
        /* 'waveletTransform3D:4' B=waveletDecomposition3D(A, w1, w2); */
        waveletDecomposition3D(A, b_w1, b_w2, B);
    } else if (direction == -1.0) {
        /* 'waveletTransform3D:5' elseif direction==-1 */
        /* 'waveletTransform3D:6' B=waveletReconstruction3D(A, w1, w2); */
        waveletReconstruction3D(A, b_w1, b_w2, B);
    } else {
        /* 'waveletTransform3D:7' else */
        /* 'waveletTransform3D:8' B=0; */
        vlen = B->size[0] * B->size[1] * B->size[2];
        B->size[0] = 1;
        B->size[1] = 1;
        B->size[2] = 1;
        emxEnsureCapacity((emxArray__common *)B, vlen, (int32_T)sizeof(creal_T));
        B->data[0].re = 0.0;
        B->data[0].im = 0.0;
    }
    emxFree_real_T(&b_w2);
    emxFree_real_T(&b_w1);
    emlrtHeapReferenceStackLeaveFcn();
}

void waveletTransform3D_api(const mxArray * const prhs[4], const mxArray *plhs[1])
{
    emxArray_creal_T *A;
    emxArray_real_T *w1;
    emxArray_real_T *w2;
    emxArray_creal_T *B;
    real_T direction;
    emlrtHeapReferenceStackEnterFcn();
    emxInit_creal_T(&A, 3, TRUE);
    emxInit_real_T(&w1, 2, TRUE);
    emxInit_real_T(&w2, 2, TRUE);
    emxInit_creal_T(&B, 3, TRUE);
    /* Marshall function inputs */
    emlrt_marshallIn(emlrtAliasP(prhs[0]), "A", A);
    c_emlrt_marshallIn(emlrtAliasP(prhs[1]), "w1", w1);
    c_emlrt_marshallIn(emlrtAliasP(prhs[2]), "w2", w2);
    direction = e_emlrt_marshallIn(emlrtAliasP(prhs[3]), "direction");
    /* Invoke the target function */
    waveletTransform3D(A, w1, w2, direction, B);
    /* Marshall function outputs */
    plhs[0] = emlrt_marshallOut(B);
    emxFree_creal_T(&B);
    emxFree_real_T(&w2);
    emxFree_real_T(&w1);
    emxFree_creal_T(&A);
    emlrtHeapReferenceStackLeaveFcn();
}

void waveletTransform3D_atexit(void)
{
    emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void waveletTransform3D_initialize(emlrtContext *context)
{
    emlrtEnterRtStack(&emlrtContextGlobal);
    emlrtFirstTime(context);
}

void waveletTransform3D_terminate(void)
{
    emlrtLeaveRtStack(&emlrtContextGlobal);
}
/* End of code generation (waveletTransform3D.c) */
