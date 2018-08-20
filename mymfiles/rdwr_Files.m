function rdwr_Files(x_dim,y_dim,n_cal)

datadir = 'C:\Documents and Settings\ashish\Desktop\MRdata\xp_complex_data\DV-XP030\2005-07-23-14-40-05\t1_mprage_sag\acq13__t1_mprage_sag';
s=['.MR..13.81.2005.07.23.15.03.16.203125.9748718.IMA';
   '.MR..13.82.2005.07.23.15.03.16.203125.9748718.IMA';
   '.MR..13.83.2005.07.23.15.03.16.203125.9748718.IMA';
   '.MR..13.84.2005.07.23.15.03.16.203125.9748718.IMA';
   '.MR..13.85.2005.07.23.15.03.16.203125.9748718.IMA';
   '.MR..13.86.2005.07.23.15.03.16.203125.9748718.IMA';
   '.MR..13.87.2005.07.23.15.03.16.203125.9748718.IMA';
   '.MR..13.88.2005.07.23.15.03.16.203125.9748718.IMA';];
thisdir = pwd;
cd(datadir);

xdim=x_dim;
ydim=y_dim;
ncal=n_cal
c=cellstr(s);
num_coils=8;
for coil=1:num_coils
   filnam=char(c(coil));
   [kData]=rdDcmRaw(filnam,xdim,ydim);
   kData=reshape(kData,xdim,ydim,1);
   data1(:,:,coil)=kData(:,:,1);   
end
whos data1
save 'data1_new' data1
tmp(:,:)=data1(:,:,7);
F2 = fftshift(fft2(tmp,xdim,ydim));
pcolor(abs(F2)); shading flat; colormap(bone);
return

cd(thisdir);

function [kData]=rdDcmRaw(filnam,xdim,ydim)
info = imfinfo(filnam,'ras');
headSize=info.FileSize-xdim*ydim*8
id=fopen(filnam,'r');
[header, count] = fread(id, headSize, 'uchar');   
[kData_1D, count] = fread(id, xdim*ydim*2, 'float');
fclose(id);
kReal=kData_1D(1:2:end);
kImag=kData_1D(2:2:end);
kData=complex(kReal, kImag);
kData=reshape(kData,xdim,ydim);
return



