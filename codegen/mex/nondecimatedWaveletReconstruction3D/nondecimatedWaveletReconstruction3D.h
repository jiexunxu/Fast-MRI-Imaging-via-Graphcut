/*
 * nondecimatedWaveletReconstruction3D.h
 *
 * Code generation for function 'nondecimatedWaveletReconstruction3D'
 *
 * C source code generated on: Thu Sep 06 12:17:49 2012
 *
 */

#ifndef __NONDECIMATEDWAVELETRECONSTRUCTION3D_H__
#define __NONDECIMATEDWAVELETRECONSTRUCTION3D_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blascompat32.h"
#include "rtwtypes.h"
#include "nondecimatedWaveletReconstruction3D_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void nondecimatedWaveletReconstruction3D(const emxArray_creal_T *A, emxArray_creal_T *B);
extern void nondecimatedWaveletReconstruction3D_api(const mxArray * const prhs[1], const mxArray *plhs[1]);
extern void nondecimatedWaveletReconstruction3D_atexit(void);
extern void nondecimatedWaveletReconstruction3D_initialize(emlrtContext *context);
extern void nondecimatedWaveletReconstruction3D_terminate(void);
#endif
/* End of code generation (nondecimatedWaveletReconstruction3D.h) */
