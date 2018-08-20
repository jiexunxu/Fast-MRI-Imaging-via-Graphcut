/*
 * waveletTransform3D_mex.c
 *
 * Code generation for function 'waveletTransform3D'
 *
 * C source code generated on: Wed Sep 05 17:20:06 2012
 *
 */

/* Include files */
#include "mex.h"
#include "waveletTransform3D.h"

/* Type Definitions */

/* Function Declarations */
static void waveletTransform3D_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "waveletTransform3D", NULL, false, NULL };

/* Function Definitions */
static void waveletTransform3D_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Temporary copy for mex outputs. */
  mxArray *outputs[1];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  /* Check for proper number of arguments. */
  if(nrhs != 4) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:WrongNumberOfInputs","4 inputs required for entry-point 'waveletTransform3D'.");
  } else if(nlhs > 1) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:TooManyOutputArguments","Too many output arguments for entry-point 'waveletTransform3D'.");
  }
  /* Module initialization. */
  waveletTransform3D_initialize(&emlrtContextGlobal);
  /* Call the function. */
  waveletTransform3D_api(prhs,(const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  waveletTransform3D_terminate();
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(waveletTransform3D_atexit);
  emlrtClearAllocCount(&emlrtContextGlobal, 0, 0, NULL);
  /* Dispatch the entry-point. */
  waveletTransform3D_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (waveletTransform3D_mex.c) */
