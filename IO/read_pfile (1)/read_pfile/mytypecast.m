% typecast creates an exact bit copy of the input as a different class
%*************************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    typecast
%  Filename:    typecast.c
%  Programmer:  James Tursa
%  Version:     1.0
%  Date:        November 6, 2007
%  Copyright:   (c) 2007 by James Tursa
%  Permission:  Permission is granted to freely distribute and use this code as long
%               as the header information is included.
% 
%  typecast is a mex function intended to mimic the MATLAB intrinsic typecast function
%  for those users with older versions of MATLAB that do not have this intrinsic.
% 
%  Building:
% 
%  >> mex -setup
%    (then follow instructions to select a C or C++ compiler of your choice)
%  >> mex typecast.c
% 
%  The usage is as follows (from The Mathworks website documentation):
% 
%  Syntax
% 
%  Y = typecast(X, type)
% 
%  Description
% 
%  Y = typecast(X, type) converts a numeric value in X to the data type specified by type.
%  Input X must be a full, noncomplex, numeric scalar or vector. The type input is a string
%  set to one of the following: 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32',
%  'uint64', 'int64', 'single', or 'double'. typecast is different from the MATLAB cast
%  function in that it does not alter the input data. typecast always returns the same
%  number of bytes in the output Y as were in the input X. For example, casting the 16-bit
%  integer 1000 to uint8 with typecast returns the full 16 bits in two 8-bit segments
%  (3 and 232) thus keeping its original value (3*256 + 232 = 1000). The cast function,
%  on the other hand, truncates the input value to 255.
% 
%  The output of typecast can be formatted differently depending on what system you use it on.
%  Some computer systems store data starting with its most significant byte (an ordering
%  called big-endian), while others start with the least significant byte (called little-endian). 
% 
%  typecast issues an error if X contains fewer values than are needed to make an output value. 
% 
%*************************************************************************************
%
function C = mytypecast(A,B)
disp('Error using typecast: You have not yet generated the typecast mex routine.');
disp('Do the following:');
disp(' ');
disp('>> mex -setup');
disp('  (Then follow instructins to select the C compiler of your choice. (e.g., lcc)');
disp('>> mex typecast.c');
disp(' ');
disp('That''s it! Now you are ready to use typecast.');
error(' ');

