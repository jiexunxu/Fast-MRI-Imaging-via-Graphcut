/*
 * nondecimatedWaveletTransform3D_mex.c
 *
 * Code generation for function 'nondecimatedWaveletTransform3D'
 *
 * C source code generated on: Thu Sep 06 12:09:03 2012
 *
 */

/* Include files */
#include "mex.h"
#include "nondecimatedWaveletTransform3D.h"

/* Type Definitions */

/* Function Declarations */
static void nondecimatedWaveletTransform3D_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "nondecimatedWaveletTransform3D", NULL, false, NULL };

/* Function Definitions */
static void nondecimatedWaveletTransform3D_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Temporary copy for mex outputs. */
  mxArray *outputs[1];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  /* Check for proper number of arguments. */
  if(nrhs != 1) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:WrongNumberOfInputs","1 input required for entry-point 'nondecimatedWaveletTransform3D'.");
  } else if(nlhs > 1) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:TooManyOutputArguments","Too many output arguments for entry-point 'nondecimatedWaveletTransform3D'.");
  }
  /* Module initialization. */
  nondecimatedWaveletTransform3D_initialize(&emlrtContextGlobal);
  /* Call the function. */
  nondecimatedWaveletTransform3D_api(prhs,(const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  nondecimatedWaveletTransform3D_terminate();
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(nondecimatedWaveletTransform3D_atexit);
  emlrtClearAllocCount(&emlrtContextGlobal, 0, 0, NULL);
  /* Dispatch the entry-point. */
  nondecimatedWaveletTransform3D_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (nondecimatedWaveletTransform3D_mex.c) */
