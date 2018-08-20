@echo off
set MATLAB=X:\PROGRA~1\MATLAB~1
set MATLAB_ARCH=win64
set MATLAB_BIN="X:\Program Files\MATLAB_R2011a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=X:\Projects\ApproximateNullVec\deployment\codegen\mex\waveletTransform3D\
set LIB_NAME=waveletTransform3D_mex
set MEX_NAME=waveletTransform3D_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for waveletTransform3D > waveletTransform3D_mex.mki
echo COMPILER=%COMPILER%>> waveletTransform3D_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> waveletTransform3D_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> waveletTransform3D_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> waveletTransform3D_mex.mki
echo LINKER=%LINKER%>> waveletTransform3D_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> waveletTransform3D_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> waveletTransform3D_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> waveletTransform3D_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> waveletTransform3D_mex.mki
echo BORLAND=%BORLAND%>> waveletTransform3D_mex.mki
echo OMPFLAGS= >> waveletTransform3D_mex.mki
echo EMC_COMPILER=msvc100free>> waveletTransform3D_mex.mki
echo EMC_CONFIG=optim>> waveletTransform3D_mex.mki
"X:\Program Files\MATLAB_R2011a\bin\win64\gmake" -B -f waveletTransform3D_mex.mk
