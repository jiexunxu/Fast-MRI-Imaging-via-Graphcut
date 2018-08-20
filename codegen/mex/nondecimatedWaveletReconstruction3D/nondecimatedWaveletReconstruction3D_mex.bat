@echo off
set MATLAB=X:\PROGRA~1\MATLAB~1
set MATLAB_ARCH=win64
set MATLAB_BIN="X:\Program Files\MATLAB_R2011a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=X:\Projects\ApproximateNullVec\deployment\codegen\mex\nondecimatedWaveletReconstruction3D\
set LIB_NAME=nondecimatedWaveletReconstruction3D_mex
set MEX_NAME=nondecimatedWaveletReconstruction3D_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for nondecimatedWaveletReconstruction3D > nondecimatedWaveletReconstruction3D_mex.mki
echo COMPILER=%COMPILER%>> nondecimatedWaveletReconstruction3D_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> nondecimatedWaveletReconstruction3D_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> nondecimatedWaveletReconstruction3D_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> nondecimatedWaveletReconstruction3D_mex.mki
echo LINKER=%LINKER%>> nondecimatedWaveletReconstruction3D_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> nondecimatedWaveletReconstruction3D_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> nondecimatedWaveletReconstruction3D_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> nondecimatedWaveletReconstruction3D_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> nondecimatedWaveletReconstruction3D_mex.mki
echo BORLAND=%BORLAND%>> nondecimatedWaveletReconstruction3D_mex.mki
echo OMPFLAGS= >> nondecimatedWaveletReconstruction3D_mex.mki
echo EMC_COMPILER=msvc100free>> nondecimatedWaveletReconstruction3D_mex.mki
echo EMC_CONFIG=optim>> nondecimatedWaveletReconstruction3D_mex.mki
"X:\Program Files\MATLAB_R2011a\bin\win64\gmake" -B -f nondecimatedWaveletReconstruction3D_mex.mk
