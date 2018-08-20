@echo off
set MATLAB=X:\PROGRA~1\MATLAB~1
set MATLAB_ARCH=win64
set MATLAB_BIN="X:\Program Files\MATLAB_R2011a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=X:\Projects\ApproximateNullVec\deployment\codegen\mex\nondecimatedWaveletTransform3D\
set LIB_NAME=nondecimatedWaveletTransform3D_mex
set MEX_NAME=nondecimatedWaveletTransform3D_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for nondecimatedWaveletTransform3D > nondecimatedWaveletTransform3D_mex.mki
echo COMPILER=%COMPILER%>> nondecimatedWaveletTransform3D_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> nondecimatedWaveletTransform3D_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> nondecimatedWaveletTransform3D_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> nondecimatedWaveletTransform3D_mex.mki
echo LINKER=%LINKER%>> nondecimatedWaveletTransform3D_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> nondecimatedWaveletTransform3D_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> nondecimatedWaveletTransform3D_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> nondecimatedWaveletTransform3D_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> nondecimatedWaveletTransform3D_mex.mki
echo BORLAND=%BORLAND%>> nondecimatedWaveletTransform3D_mex.mki
echo OMPFLAGS= >> nondecimatedWaveletTransform3D_mex.mki
echo EMC_COMPILER=msvc100free>> nondecimatedWaveletTransform3D_mex.mki
echo EMC_CONFIG=optim>> nondecimatedWaveletTransform3D_mex.mki
"X:\Program Files\MATLAB_R2011a\bin\win64\gmake" -B -f nondecimatedWaveletTransform3D_mex.mk
