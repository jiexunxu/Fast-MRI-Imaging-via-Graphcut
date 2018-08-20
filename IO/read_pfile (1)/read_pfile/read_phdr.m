function phdr = read_phdr( pfilename )
%GEHDR_READER Summary of this function goes here
%   Detailed explanation goes here
esever=getESE(pfilename);
gehdr_sizeof=str2func(['gehdr_sizeof' esever]);
gehdr_reader=str2func(['gehdr_reader' esever]);
headersize=gehdr_sizeof('POOL_HEADER');
fid=fopen(pfilename,'rb');chr=fread(fid,headersize,'*uchar');fclose(fid);
phdr=gehdr_reader(chr, 1, 'POOL_HEADER', 1);
end

