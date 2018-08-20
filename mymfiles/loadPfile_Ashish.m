function zftdata=loadPfile_Ashish(pfile, xres, yres, nzpe, nc, headersize, type, mcskip, PA_offset)
try
    if nargin==9
        zftdata=tmp_loadPfile_Ashish(pfile, nc, nzpe, xres, yres, headersize, type, mcskip, PA_offset);
    elseif nargin==8
        zftdata=tmp_loadPfile_Ashish(pfile, nc, nzpe, xres, yres, headersize, type, mcskip);
    elseif nargin==7
        zftdata=tmp_loadPfile_Ashish(pfile, nc, nzpe, xres, yres, headersize, type);
    elseif nargin==6
        zftdata=tmp_loadPfile_Ashish(pfile, nc, nzpe, xres, yres, headersize);
    elseif nargin==5
        zftdata=tmp_loadPfile_Ashish(pfile, nc, nzpe, xres, yres);
    elseif nargin==4
        zftdata=tmp_loadPfile_Ashish(pfile, 1, nzpe, xres, yres);
    elseif nargin==3
        zftdata=tmp_loadPfile_Ashish(pfile, 1, 1, xres, yres);
    end
catch
    try
        zftdata=tmp_loadPfile_Ashish(pfile, xres, yres, nzpe, nc, 66072);
    catch
        try
            zftdata=tmp_loadPfile_Ashish(pfile, xres, yres, nzpe, nc, 61464);
        catch
            try
                zftdata=tmp_loadPfile_Ashish(pfile, xres, yres, nzpe, nc, 39984);
            catch
                zftdata=tmp_loadPfile_Ashish(pfile, xres, yres, nzpe, nc, 145908);
            end
        end
    end
end
function zftdata=tmp_loadPfile_Ashish(pfile, xres, yres, nzpe, nc, headersize, type, mcskip, PA_offset)
%% By Ashish Raj - combines all knwon Pfile reader code into one - tries each one in succession and stops if read successful
%zftdata = zeros(xres, yres, nzpe, nc, 'single');
try
    if nargin==9
        zftdata=loadPfile_new(pfile, nc, nzpe, xres, yres, headersize, type, mcskip, PA_offset);
    elseif nargin==8
        zftdata=loadPfile_new(pfile, nc, nzpe, xres, yres, headersize, type, mcskip);
    elseif nargin==7
        zftdata=loadPfile_new(pfile, nc, nzpe, xres, yres, headersize, type);
    elseif nargin==6
        zftdata=loadPfile_new(pfile, nc, nzpe, xres, yres, headersize);
    elseif nargin==5
        zftdata=loadPfile_new(pfile, nc, nzpe, xres, yres);
    elseif nargin==4
        zftdata=loadPfile_new(pfile, 1, nzpe, xres, yres);
    elseif nargin==3
        zftdata=loadPfile_new(pfile, 1, 1, xres, yres);
    else
        disp('not enough input parameters!');
    end  
catch
    disp('could not read using loadPfile_new.m... Now trying loadPfile_Ramin.m ');
    %clear zftdata;
    try
        if nargin>=9
            zftdata=loadPfile_Ramin(pfile, nc, nzpe, xres, yres, headersize, type, mcskip, PA_offset);
        elseif nargin>=8
            zftdata=loadPfile_Ramin(pfile, nc, nzpe, xres, yres, headersize, type, mcskip);
        elseif nargin>=7
            zftdata=loadPfile_Ramin(pfile, nc, nzpe, xres, yres, headersize, type);
        elseif nargin>=6
            zftdata=loadPfile_Ramin(pfile, nc, nzpe, xres, yres, headersize);
        elseif nargin>=5
            zftdata=loadPfile_Ramin(pfile, nc, nzpe, xres, yres);
        elseif nargin>=4
            zftdata=loadPfile_Ramin(pfile, 1, nzpe, xres, yres);
        elseif nargin>=3
            zftdata=loadPfile_Ramin(pfile, 1, 1, xres, yres);
        else
            disp('not enough input parameters!');
        end
    catch
        disp('could not read using loadPfile_Ramin.m... Now trying loadPfile_Thanh.m');
        %clear zftdata;
        try
            if nargin>=6
                    zftdata=loadPfile_Thanh(pfile,xres,yres,nzpe,nc, headersize);
            elseif nargin>=5
                    zftdata=loadPfile_Thanh(pfile,xres,yres,nzpe,nc);
            elseif nargin>=4
                    zftdata=loadPfile_Thanh(pfile,xres,yres,nzpe,1);
            elseif nargin>=3
                    zftdata=loadPfile_Thanh(pfile,xres,yres,1,1);
            else
                disp('not enough input parameters!');
            end
        catch
            disp('could not read using loadPfile_Thanh.m... Now trying ToAshish_readPfile.m');
            %clear zftdata;
            try
                if nargin>=6
                        zftdata=ToAshish_readPfile(pfile, xres, yres, nzpe, nc, headersize);
                elseif nargin>=5
                        zftdata=ToAshish_readPfile(pfile, xres, yres, nzpe, nc);
                elseif nargin>=4
                        zftdata=ToAshish_readPfile(pfile, xres, yres, nzpe, 1);
                elseif nargin>=3
                        zftdata=ToAshish_readPfile(pfile, xres, yres, 1, 1);
                else
                    disp('not enough input parameters!');
                end       
            catch
                disp('Tried several Pfile readers, but failed to read raw data!');
            end
        end
    end
end
read_volume_size = size(zftdata),
zftdata = squeeze(zftdata);

%% Internal functions - is simply cut and pasted from elsewhere - but this is intended as the most current version

function zftdata=loadPfile_new(pfile, nc, nzpe, xres, yres, headersize, type, mcskip, PA_offset)
    %input -- pfile: p-file name and directory
    %      -- nc: the number of coils
    %      -- nzpe: number of slices or z phase encoding
    %      -- xres: the frequency encoding number, for spiral it's the leaf length (ndat)
    %      -- yres: the phase encoding number, for spiral it's the number of leaves (nl)
    %      -- type: psd type: 0-conventional, 1-spiral
    %      -- baseline: the number of bytes for baseline
    %      -- headersize: the number of bytes for header
    %output-- zftdata: the rawdata after fourier transform in z direction


    if nargin < 9, PA_offset = 20; end
    if nargin < 8, mcskip = 0; end
    if nargin < 7, type = 0; end
    
%     % For Version 7 Pfiles
%     if nargin < 6  || isempty(headersize), headersize = 39984; end
    % For Version 9 Pfiles
    % if nargin < 6  || isempty(headersize), headersize = 61464; end
    if nargin < 6  || isempty(headersize), headersize = 66072; end
    
    if nargin < 5, yres=256; end
    if nargin < 4, xres=256; end
    if nargin < 3, nzpe = 1; end
    if nargin < 2, nc = 4; end

    %nc, nzpe, xres, yres, PA_offset, headersize,

% below endian verification does not seem to  work - replacing by little
% endian
%     fid = fopen(pfile, 'r', 'ieee-le');
%     leversion = fread(fid,1,'float32');
%     fclose(fid);
% 
%     fid = fopen(pfile, 'r', 'ieee-be');
%     beversion = fread(fid,1,'float32');
%     fclose(fid);
%     %leversion, beversion,
%     endian_str = 'ieee-le';
% 
%     %fseek(fid, 0, 'bof');           % Not closed because used below
%     if (beversion == 7)
%          disp('Version 7 Pfile');
%          endian_str = 'ieee-be';
%     end
% 
%     if (leversion == 9)
%          disp('Version 9 Pfile');
%          endian_str = 'ieee-le';
%     end

    endian_str = 'ieee-le';
    fid = fopen(pfile, 'r', endian_str);

    if (type==0)
        status = fseek(fid, headersize, 'bof');
        [raw, c]=fread(fid, 'int16');
        status=fclose(fid);
        Re(1:c/2)=raw(1:2:c-1);
        Im(1:c/2)=raw(2:2:c);
        zftdata=complex(Re,Im);
        clear Re Im raw;
        zftdata=reshape(zftdata, xres, yres+1, nzpe, nc);
        zftdata=zftdata(1:xres, 2:yres+1, 1:nzpe, 1:nc);
        %zftdata(1:xres, 1:yres, 1:nzpe, 1:nc)=zftdata(1:xres, 2:yres+1, 1:nzpe, 1:nc);
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
            prd3d=complex(Re,Im);
            clear Re Im raw;
            prd3d=reshape(prd3d, ndat, nl, nzpe);

            %fft in kz direction
            zft3d=zeros(ndat,nl,nzpe);
            zft3d=fftshift(fft(prd3d, nzpe, 3),3);
            zft3d(ndat-nskip:ndat,:,:)=0;
            zftdata(:,:,:,coilnum)=zft3d;
        end
        status=fclose(fid);
    end
    zftdata = single(zftdata);
end

function zftdata=loadPfile_Ramin(pfile, nc, nzpe, xres, yres, headersize, type, mcskip, PA_offset)

    %input -- pfile: p-file name and directory
    %      -- nc: the number of coils
    %      -- nzpe: number of slices or z phase encoding
    %      -- xres: the frequency encoding number, for spiral it's the leaf
    %      length (ndat)
    %      -- yres: the phase encoding number, for spiral it's the number of
    %      leaves (nl)
    %      -- type: psd type: 0-conventional, 1-spiral
    %      -- baseline: the number of bytes for baseline
    %output-- zftdata: the rawdata after fourier transform in z direction

    if nargin < 9, PA_offset = 20; end
    if nargin < 8, mcskip = 0; end
    if nargin < 7, type = 0; end
    if nargin < 5, yres=256; end
    if nargin < 4, xres=256; end
    if nargin < 3, nzpe = 1; end
    if nargin < 2, nc = 4; end

    fid = fopen(pfile, 'r', 'ieee-le');
    leversion = fread(fid,1,'float32');
    fclose(fid);
    fid = fopen(pfile, 'r', 'ieee-be');
    beversion = fread(fid,1,'float32');
    %leversion, beversion,
    fclose(fid);
    if (leversion == 9)
         disp('Version 9 Pfile');
         fmt = 'int32';
         ordr = 'ieee-le';
         if exist('headersize')~=1 || isempty(headersize),
             headersize = 61464;
         end
    elseif (beversion == 7)
         disp('Version 7 Pfile');
         fmt = 'int16';
         ordr = 'ieee-be';
          if exist('headersize')~=1  || isempty(headersize),
             headersize = 39984;
          end
    else
         disp('Unknown Pfile type');
    end
    fid = fopen(pfile, 'r', ordr);
    if (type==0)
        status = fseek(fid, headersize, 'bof');
        [raw, c]=fread(fid, fmt);
        status=fclose(fid);
        Re(1:c/2)=raw(1:2:c-1);
        Im(1:c/2)=raw(2:2:c);
        raw_c=complex(Re,Im);
        clear Re Im raw;
        zftdata=single(reshape(raw_c, xres, yres+1, nzpe, nc));
        clear raw_c;
        zftdata = zftdata(1:xres, 2:yres+1, 1:nzpe, 1:nc);
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
            clear Re Im raw;
            prd3d=zeros(nzpe,nl,ndat);
            prd3d=reshape(raw_c, ndat, nl, nzpe);
            clear raw_c;
            %fft in kz direction
            zft3d=zeros(ndat,nl,nzpe);
            zft3d=fftshift(fft(prd3d, nzpe, 3),3);
            zft3d(ndat-nskip:ndat,:,:)=0;
            zftdata(:,:,:,coilnum)=single(zft3d);
        end
        status=fclose(fid);
    end
    zftdata = single(zftdata);
end

function imf = loadPfile_Thanh(filename,Nx,Ny,Nz,numrecv, headersize)
    % was found successful in reading Thanh's 2D cardiac cine noise data

    if nargin<6 || isempty(headersize),
        offset = 145908;
    else
        offset =  headersize;
    end
    fid = fopen(filename,'r','b');
    stat = fseek(fid,offset,'bof');         % skips header
    if stat ~= 0
    msg = ferror(fid);
    error(msg);
    end
    [data,count] = fread(fid,[Nx*2,inf],'int16');
    fclose(fid);
    imf = zeros(Nx, Ny, Nz, numrecv, 'single');
    for m = 0:numrecv-1
    for n = 0:Nz-1
        tmp = data(:,2+n*(Ny+1)+m*Nz*(Ny+1):Ny+1+n*(Ny+1)+m*Nz*(Ny+1)); % skips baseline
        tmp = tmp(1:2:size(tmp,1),:) + sqrt(-1)*tmp(2:2:size(tmp,1),:);
        imf(:, :, n+1, m+1) = single(tmp);
    end
    end
    imf = single(imf);
end

function rawmatt = ToAshish_readPfile(pfile, freq, phase, slices, ncoils, headersize)
    if nargin<6 || isempty(headersize),
        HeaderSize = 66072;
    else
        HeaderSize =  headersize;
    end
    % READING PFILE, MODULATION AND SIZE CORRECTION
    fid = fopen(pfile, 'r');
    fseek(fid, HeaderSize, -1);
    raw = fread(fid, inf, 'int32');
    fclose (fid);
    rawc = raw(1:2:end) + j.*raw(2:2:end);
    clear raw;
    rawmatt = single(reshape(rawc, freq, phase + 1, slices, ncoils));    %fixed 8 for 8 coils
    clear rawc;
    rawmatt = rawmatt(:,2:end,:,:);
    rawmatt = single(rawmatt);
end

end  % of main function

    
end