function [ data hdr ] = read_pfile( pfile )
%PFILE_READ Summary of this function goes here
%   Detailed explanation goes here
hdr=read_phdr(pfile);
offset=size_phdr(pfile);
xres=double(hdr.rdb_hdr_rec.rdb_hdr_da_xres);
yres=double(hdr.rdb_hdr_rec.rdb_hdr_da_yres);
zres=double(hdr.rdb_hdr_rec.rdb_hdr_nslices);
tres=double(hdr.rdb_hdr_rec.rdb_hdr_dab.stop_rcv(1))-...
    double(hdr.rdb_hdr_rec.rdb_hdr_dab.start_rcv(1))+1;
fid=fopen(pfile);fseek(fid,offset,'bof');
data=fread(fid,2*xres*yres*zres*tres,'short');
data=complex(data(1:2:end),data(2:2:end));
data=reshape(data, [xres yres zres tres]);
data=data(:,2:end,:,:);
fclose(fid);
end

