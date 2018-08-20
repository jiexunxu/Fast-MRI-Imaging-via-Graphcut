/*
 * waveletTransform3D.h
 *
 * Code generation for function 'waveletTransform3D'
 *
 * C source code generated on: Wed Sep 05 17:20:06 2012
 *
 */

#ifndef __WAVELETTRANSFORM3D_H__
#define __WAVELETTRANSFORM3D_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blascompat32.h"
#include "rtwtypes.h"
#include "waveletTransform3D_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void waveletTransform3D(const emxArray_creal_T *A, const emxArray_real_T *w1, const emxArray_real_T *w2, real_T direction, emxArray_creal_T *B);
extern void waveletTransform3D_api(const mxArray * const prhs[4], const mxArray *plhs[1]);
extern void waveletTransform3D_atexit(void);
extern void waveletTransform3D_initialize(emlrtContext *context);
extern void waveletTransform3D_terminate(void);
#endif
/* End of code generation (waveletTransform3D.h) */
