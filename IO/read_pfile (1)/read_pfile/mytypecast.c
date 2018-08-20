/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    typecast
 * Filename:    typecast.c
 * Programmer:  James Tursa
 * Version:     1.1
 * Date:        November 6, 2007
 * Copyright:   (c) 2007 by James Tursa
 * Permission:  Permission is granted to freely distribute and use this code as long
 *              as the header information is included.
 *
 * typecast is a mex function intended to mimic the MATLAB intrinsic typecast function
 * for those users with older versions of MATLAB that do not have this intrinsic.
 *
 * Building:
 *
 * >> mex -setup
 *   (then follow instructions to select a C or C++ compiler of your choice)
 * >> mex typecast.c
 *
 * The usage is as follows (from The Mathworks website documentation):
 *
 * Syntax
 *
 * Y = typecast(X, type)
 *
 * Description
 *
 * Y = typecast(X, type) converts a numeric value in X to the data type specified by type.
 * Input X must be a full, noncomplex, numeric scalar or vector. The type input is a string
 * set to one of the following: 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32',
 * 'uint64', 'int64', 'single', or 'double'. typecast is different from the MATLAB cast
 * function in that it does not alter the input data. typecast always returns the same
 * number of bytes in the output Y as were in the input X. For example, casting the 16-bit
 * integer 1000 to uint8 with typecast returns the full 16 bits in two 8-bit segments
 * (3 and 232) thus keeping its original value (3*256 + 232 = 1000). The cast function,
 * on the other hand, truncates the input value to 255.
 *
 * The output of typecast can be formatted differently depending on what system you use it on.
 * Some computer systems store data starting with its most significant byte (an ordering
 * called big-endian), while others start with the least significant byte (called little-endian). 
 *
 * typecast issues an error if X contains fewer values than are needed to make an output value. 
 *
 */

#include "mex.h"
#include "matrix.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    int inbytes, outbytes, m, n, k, d, numel;
    size_t inbytes_t, numel_t, nbytes_t;
    char *outstring;
    mxClassID outclass;

/* Check input arguments */
    
   if( nrhs > 2 )
       mexErrMsgTxt("Too many input arguments.\n");
   if( nrhs < 2 )
       mexErrMsgTxt("Not enough input arguments.\n");
   if( !(mxIsNumeric(prhs[0])||mxIsChar(prhs[0])) || mxIsComplex(prhs[0]) || mxIsEmpty(prhs[0]) )
       mexErrMsgTxt("The first input argument must be a full, non-complex numeric value.\n");
   if( !mxIsChar(prhs[1]) )
       mexErrMsgTxt("The second input argument must be a character array.\n");
   if( mxGetNumberOfDimensions(prhs[0]) != 2 )
       mexErrMsgTxt("The first input argument must be a vector.\n");
   m = mxGetM(prhs[0]);
   n = mxGetN(prhs[0]);
   if( m != 1 && n != 1 )
       mexErrMsgTxt("The first input argument must be a vector.\n");

/* Check first input argument for byte length */

   switch( mxGetClassID(prhs[0]) )
   {
       case mxINT8_CLASS:
       case mxUINT8_CLASS:
       case mxCHAR_CLASS:
           inbytes = 1;
           break;
       
       case mxINT16_CLASS:
       case mxUINT16_CLASS:
           inbytes = 2;
           break;
    
       case mxINT32_CLASS:
       case mxUINT32_CLASS:
       case mxSINGLE_CLASS:
           inbytes = 4;
           break;

       case mxINT64_CLASS:
       case mxUINT64_CLASS:
       case mxDOUBLE_CLASS:
           inbytes = 8;
           break;
       
       default:
           mexErrMsgTxt("The first input argument is not a supported type.\n");
       break;
   }
   
/* Check second input argument for desired output type */
   
   outstring = mxArrayToString(prhs[1]);
   
   if( strcmp(outstring,"int8") == 0 )
   {
       outclass = mxINT8_CLASS;
       outbytes = 1;
   }
   else if( strcmp(outstring,"uint8") == 0 )
   {
       outclass = mxUINT8_CLASS;
       outbytes = 1;
   }
   else if( strcmp(outstring,"int16") == 0 )
   {
       outclass = mxINT16_CLASS;
       outbytes = 2;
   }
   else if( strcmp(outstring,"uint16") == 0 )
   {
       outclass = mxUINT16_CLASS;
       outbytes = 2;
   }
   else if( strcmp(outstring,"int32") == 0 )
   {
       outclass = mxINT32_CLASS;
       outbytes = 4;
   }
   else if( strcmp(outstring,"uint32") == 0 )
   {
       outclass = mxUINT32_CLASS;
       outbytes = 4;
   }
   else if( strcmp(outstring,"int64") == 0 )
   {
       outclass = mxINT64_CLASS;
       outbytes = 8;
   }
   else if( strcmp(outstring,"uint64") == 0 )
   {
       outclass = mxUINT64_CLASS;
       outbytes = 8;
   }
   else if(( strcmp(outstring,"double") == 0 ) || ( strcmp(outstring,"float64") == 0 ))
   {
       outclass = mxDOUBLE_CLASS;
       outbytes = 8;
   }
   else if(( strcmp(outstring,"single") == 0 ) || ( strcmp(outstring,"float32") == 0 ) )
   {
       outclass = mxSINGLE_CLASS;
       outbytes = 4;
   }
   else if(( strcmp(outstring,"uchar") == 0 ) || ( strcmp(outstring,"schar") == 0 ))
   {
       outclass = mxCHAR_CLASS;
       outbytes = 1;
   }
   else
   {
       mxFree(outstring);
       mexErrMsgTxt("Unsupported class.\n");
   }
   mxFree(outstring);

/* Number of elements to copy*/
   
   numel = m > n ? m : n;
   numel_t = numel;
   
/* Calculate new row or column value, check for overflow & matching lengths*/
   
   if( inbytes > outbytes )
   {
       k = numel * (inbytes / outbytes);
       if( k < numel )
           mexErrMsgTxt("Array too large.\n");
       numel = k;
   }
   else
   {
       d = (outbytes / inbytes);
       k = numel / d;
       if( k * d != numel )
           mexErrMsgTxt("Too few input values to make output type.\n");
       numel = k;
   }

/* Calculate the number of bytes to copy, check for overflow*/
   
   inbytes_t = inbytes;
   nbytes_t = inbytes_t * numel_t;
   if( nbytes_t < numel_t && numel_t != 0 )
       mexErrMsgTxt("Array too large.\n");
   
/* Create the output array*/
   
   if( m > n )
       m = numel;
   else
       n = numel;
   plhs[0] = mxCreateNumericMatrix(m, n, outclass, mxREAL);

/* Copy the data*/
   
   memcpy(mxGetData(plhs[0]), mxGetData(prhs[0]), nbytes_t);

}
