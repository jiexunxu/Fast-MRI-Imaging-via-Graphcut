/*
 * nondecimatedWaveletTransform3D.h
 *
 * Code generation for function 'nondecimatedWaveletTransform3D'
 *
 * C source code generated on: Thu Sep 06 12:09:03 2012
 *
 */

#ifndef __NONDECIMATEDWAVELETTRANSFORM3D_H__
#define __NONDECIMATEDWAVELETTRANSFORM3D_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blascompat32.h"
#include "rtwtypes.h"
#include "nondecimatedWaveletTransform3D_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void nondecimatedWaveletTransform3D(const emxArray_creal_T *A, emxArray_creal_T *B);
extern void nondecimatedWaveletTransform3D_api(const mxArray * const prhs[1], const mxArray *plhs[1]);
extern void nondecimatedWaveletTransform3D_atexit(void);
extern void nondecimatedWaveletTransform3D_initialize(emlrtContext *context);
extern void nondecimatedWaveletTransform3D_terminate(void);
#endif
/* End of code generation (nondecimatedWaveletTransform3D.h) */
