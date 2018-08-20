function headersize = size_phdr( pfilename )
%GEHDR_READER Summary of this function goes here
%   Detailed explanation goes here
esever=getESE(pfilename);
gehdr_sizeof=str2func(['gehdr_sizeof' esever]);
headersize=gehdr_sizeof('POOL_HEADER');
end

