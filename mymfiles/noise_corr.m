%function noise_corr()
% analyzes the noise (no load) in multichannel parallel imaging experiments
% input: no-load noise scans, output: noise correlation matrix 

% pfile = 'L:\MR_Data\thanh_noise_data\P55808.7';
% pfile = 'L:\MR_Data\thanh_noise_data\P56832.7';
% pfile = 'L:\MR_Data\thanh_noise_data\P58368.7';
 pfile = 'L:\MR_Data\thanh_noise_data\P59392.7';

xres = 256; yres = 256; nzpe = 16; nc = 8;
zftdata=loadPfile_Ashish(pfile, xres, yres, nzpe, nc);
q = fftshift(fft2(zftdata(:,:,1,1)));
qq = abs(q);
figure; colormap(gray); imshow(qq/max(qq(:)));

[nrows, ncols, nslices, ncoils] = size(zftdata);

img = zeros(size(zftdata), 'single');
for i = 1:nslices
    for j = 1:ncoils
        img(:,:,i,j) = single(fftshift(ifft2(zftdata(:,:,i,j))));
    end
end

% spatial correlation : noise color
for k = 1:ncoils
    nvec{k} = reshape(squeeze(img(:,:,:,k)), [nrows*ncols, nslices]);
    %nvec{k} = reshape(squeeze(img(:,100,:,k)), [nrows, nslices]);
    %Csp{k} = (nvec{k}/sqrt(nslices))*(nvec{k}/sqrt(nslices))';
    spect{k} = pwelch(nvec{k}(:));
end
figure; 
for k=1:ncoils
    subplot(2,4,k); plot(spect{k});
end

% spatial correlation : noise variance

nhood = ones(41,41);
figure; 
for k=1:ncoils
    q = squeeze(abs(img(:,:,6,k)));
    mxq = max(max(q));
    q(abs(q) == mxq) = 0;  % removing DC artifact
    qq = stdfilt(q, nhood);
    subplot(2,4,k); image(qq*6); colormap(gray);
end


% cross coil correlation
% for k = 1:ncoils
%     %Csp{k} = (nvec{k}/sqrt(nslices))*(nvec{k}/sqrt(nslices))';
%     spect{k} = pwelch(nvec{k}(:));
% end

C = cov(reshape(img, [nrows*ncols*nslices, ncoils]));
figure; imagesc(abs(C)); colormap(gray); title('Cross-coil covariance matrix');
% for k=1:ncoils
%     subplot(2,4,k); plot(spect{k});
% end

