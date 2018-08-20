function zftdata = read_Pfile_with_hdr(pfile, bytesperreal, endian_str)
% uses header to determine image size and offset etc
% also does xchop etc 

if nargin<2
    bytesperreal = 2; % or 2 or 8
end
if nargin<3
    endian_str = 'ieee-le';  % or 'ieee-be';
end

if bytesperreal == 2
    fmt = 'int16';
elseif bytesperreal == 4
    fmt = 'int32';
elseif bytesperreal == 8  
    fmt = 'int64';
end
fid = fopen(pfile, 'r', endian_str);
hdr = read_gehdr(fid);
offset = hdr.rdb.off_data;
if offset==0
    offset=66072;
end
nslices = hdr.rdb.nslices;
yres = hdr.rdb.nframes;
xres = hdr.rdb.frame_size;

% determine #coils
d = dir(pfile);
numvox = (d.bytes - offset)/bytesperreal/2;
noofcoils = numvox/xres/(yres+1)/nslices;
if mod(noofcoils, 1) ~=0
    disp(sprintf('number of coils (%0.2f) was not found to be integer - check!', noofcoils));
end

status = fseek(fid, offset, 'bof');
[raw, c]=fread(fid, fmt);
status=fclose(fid);
Re(1:c/2)=raw(1:2:c-1);
Im(1:c/2)=raw(2:2:c);
raw_c=complex(Re,Im);
clear Re Im raw;
zftdata=single(reshape(raw_c, xres, yres+1, nslices, noofcoils));
clear raw_c;
zftdata = squeeze(zftdata(:, 2:yres+1, :, :));

if nslices>1
    for slice = 1:nslices    
        for i = 1:noofcoils
            zftdata(:,:,slice,i) = doYchop(doXchop(zftdata(:,:,slice,i)));
        end
    end
else
    for i = 1:noofcoils
        zftdata(:,:,i) = doYchop(doXchop(zftdata(:,:,i)));
    end
end



