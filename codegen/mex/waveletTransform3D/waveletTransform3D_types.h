/*
 * waveletTransform3D_types.h
 *
 * Code generation for function 'waveletTransform3D'
 *
 * C source code generated on: Wed Sep 05 17:20:05 2012
 *
 */

#ifndef __WAVELETTRANSFORM3D_TYPES_H__
#define __WAVELETTRANSFORM3D_TYPES_H__

/* Type Definitions */
#ifndef struct_emxArray_creal_T
#define struct_emxArray_creal_T
typedef struct emxArray_creal_T
{
    creal_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_creal_T;
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

#endif
/* End of code generation (waveletTransform3D_types.h) */
