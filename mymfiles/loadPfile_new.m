function zftdata=loadPfile_new(pfile, nc, nzpe, xres, yres, type, mcskip, PA_offset, headersize)

%input -- pfile: p-file name and directory
%      -- nc: the number of coils
%      -- nzpe: number of slices or z phase encoding
%      -- xres: the frequency encoding number, for spiral it's the leaf length (ndat)
%      -- yres: the phase encoding number, for spiral it's the number of leaves (nl)
%      -- type: psd type: 0-conventional, 1-spiral
%      -- baseline: the number of bytes for baseline
%      -- headersize: the number of bytes for header
%output-- zftdata: the rawdata after fourier transform in z direction

% For Version 7 Pfiles
if nargin < 9, headersize = 39984; end
% For Version 9 Pfiles
% if nargin < 9, headersize = 61464; end
if nargin < 9, headersize = 66072; end

if nargin < 8, PA_offset = 20; end
if nargin < 7, mcskip = 0; end
if nargin < 6, type = 0; end
if nargin < 5, yres=256; end
if nargin < 4, xres=256; end
if nargin < 3, nzpe = 1; end
if nargin < 2, nc = 4; end
 
nc, nzpe, xres, yres, PA_offset, headersize,

fid = fopen(pfile, 'r', 'ieee-le');
leversion = fread(fid,1,'float32');
fclose(fid);

fid = fopen(pfile, 'r', 'ieee-be');
beversion = fread(fid,1,'float32');
fclose(fid);
%leversion, beversion,
endian_str = 'ieee-be';

%fseek(fid, 0, 'bof');           % Not closed because used below
if (beversion == 7)
     disp('Version 7 Pfile');
     endian_str = 'ieee-be';
end

if (leversion == 9)
     disp('Version 9 Pfile');
     endian_str = 'ieee-le';
end
fid = fopen(pfile, 'r', endian_str);
    
if (type==0)
    status = fseek(fid, headersize, 'bof');
    [raw, c]=fread(fid, 'int16');
    status=fclose(fid);
    Re(1:c/2)=raw(1:2:c-1);
    Im(1:c/2)=raw(2:2:c);
    raw_c=complex(Re,Im);
    %size(raw_c), xres, yres, nzpe, nc,
    data=reshape(raw_c, xres, yres+1, nzpe, nc);
    zftdata(1:xres, 1:yres, 1:nzpe, 1:nc)=data(1:xres, 2:yres+1, 1:nzpe, 1:nc);
    zftdata = squeeze(zftdata);
else 
    baseline =  nc*xres;
    nskip = 50;
    nl=yres;
    ndat=xres;
    for coilnum = 1:nc
        fseek(fid, headersize+baseline+PA_offset+mcskip*(coilnum-1), -1);
        [raw, c]=fread(fid, nl*nzpe*ndat*2, 'int16');
        Re(1:c/2)=raw(1:2:c-1);
        Im(1:c/2)=raw(2:2:c);
        raw_c=complex(Re,Im);
        
        prd3d=zeros(nzpe,nl,ndat);
        prd3d=reshape(raw_c, ndat, nl, nzpe);
        
        %fft in kz direction
        zft3d=zeros(ndat,nl,nzpe);
        zft3d=fftshift(fft(prd3d, nzpe, 3),3);
        zft3d(ndat-nskip:ndat,:,:)=0;
        zftdata(:,:,:,coilnum)=zft3d;
    end
	status=fclose(fid);
end
