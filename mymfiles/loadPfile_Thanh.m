function imf = loadPfile_Thanh(filename,Nx,Ny,Nz,numrecv)
% was found successful in reading Thanh's 2D cardiac cine noise data

offset = 145908;
fid = fopen(filename,'r','b');

stat = fseek(fid,offset,'bof');         % skips header
if stat ~= 0
msg = ferror(fid);
error(msg);
end
[data,count] = fread(fid,[Nx*2,inf],'int16');
fclose(fid);

imf = zeros(Nx, Ny, Nz, numrecv);
for m = 0:numrecv-1
for n = 0:Nz-1
    tmp = data(:,2+n*(Ny+1)+m*Nz*(Ny+1):Ny+1+n*(Ny+1)+m*Nz*(Ny+1)); % skips baseline
    tmp = tmp(1:2:size(tmp,1),:) + i*tmp(2:2:size(tmp,1),:);
    imf(:, :, n+1, m+1) = tmp;
end
end
read_volume_size = size(imf),
