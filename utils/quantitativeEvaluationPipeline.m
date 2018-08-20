function quantitativeEvaluationPipeline(dataFileName, Gs, noises)
    data0=readAndSimulateUCSFData(dataFileName, []);
    data0=data0(:, :, 110:112);
    m=size(data0, 1);n=size(data0, 2);    
    [data, sensit]=readAndSimulateUCSFData(data0, true(m, n));
    xRef=computeSENSE(data, sensit, true(m, n), 0.25, 'noTV');    
    voxelMask=abs(xRef)>0.7*mean(abs(xRef(:)));
    voxelMask(:, :, 2)=bwmorph(voxelMask(:, :, 2), 'erode', 2);
    voxelMask(:, :, 2)=bwmorph(voxelMask(:, :, 2), 'dilate', 5);
    voxelMask(:, :, 2)=bwmorph(voxelMask(:, :, 2), 'fill');
    voxelMask(:, :, 2)=bwmorph(voxelMask(:, :, 2), 'erode', 3);
    xRef=xRef(:, :, 2).*voxelMask(:, :, 2);
    xRef=uint8(256*abs(xRef)/max(abs(xRef(:))));
    for i=1:length(Gs)
        G=Gs{i};
        [data, sensit]=readAndSimulateUCSFData(data0(:, :, 110:112), true(m, n));
        xSENSE=computeSENSE(data, sensit, G, 0.25, 'noTV');
        x0=estimateInitialSolution(data, sensit, G, 1e-9, 0.01, 100);
        x=subspaceGraphcutPipeline(data(:, :, 2, :), sensit(:, :, 2, :), G, x0(:, :, 2), voxelMask(:, :, 2));        
        xEval{i, 1}=xSENSE(:, :, 2).*voxelMask(:, :, 2);xEval{i, 2}=x0(:, :, 2).*voxelMask(:, :, 2);xEval{i, 3}=x.*voxelMask(:, :, 2);
        xEval{i, 1}=uint8(256*abs(xEval{i, 1})/max(abs(xEval{i, 1}(:))));
        xEval{i, 2}=uint8(256*abs(xEval{i, 2})/max(abs(xEval{i, 2}(:))));
        xEval{i, 3}=uint8(256*abs(xEval{i, 3})/max(abs(xEval{i, 3}(:))));        
    end
    
    for i=1:length(noises)
        [data, sensit]=readAndSimulateUCSFData(data0+noises(i)*randn(size(data0)), G);
        xSENSE=computeSENSE(data, sensit, G, 0.25, 'noTV');
        x0=estimateInitialSolution(data, sensit, G, 1e-9, 0.01, 100);
        x=subspaceGraphcutPipeline(data(:, :, 2, :), sensit(:, :, 2, :), G, x0(:, :, 2), voxelMask(:, :, 2));
        [x, xCS]=subspaceGraphcutPipeline(data, sensit, G, []);
        xResult{i+1, 1}=xSENSE;xResult{i+1, 2}=xCS;xResult{i+1, 3}=x;
        xEval2{i, 1}=xSENSE(:, :, 2).*voxelMask(:, :, 2);xEval{i, 2}=x0(:, :, 2).*voxelMask(:, :, 2);xEval{i, 3}=x.*voxelMask(:, :, 2);
        xEval2{i, 1}=uint8(256*abs(xEval2{i, 1})/max(abs(xEval2{i, 1}(:))));
        xEval2{i, 2}=uint8(256*abs(xEval2{i, 2})/max(abs(xEval2{i, 2}(:))));
        xEval2{i, 3}=uint8(256*abs(xEval2{i, 3})/max(abs(xEval2{i, 3}(:))));  
    end
    
    accScores=scoring(xEval, xRef);
    noiseScores=scoring(xEval2, xRef);
end

function scores=scoring(xEval, xRef)
    scores=zeros(size(xEval, 1), 3, 5);
    for i=1:size(xEval, 1)
        for j=1:3
            x=xEval{i, j};
            PSNRScore=psnr(xRef, x);
            [SSIMScore, ~, LSR_map]=my_ssim(x, x0, [0.01 0.03 0.015], ones(11), max(x0(:)), [1.5 0.6 0.6]);
            SSIM_L=mean2(LSR_map(:, :, 1));SSIM_S=mean2(LSR_map(:, :, 2));SSIM_R=mean2(LSR_map(:, :, 3));
            scores(i, j, 1)=PSNRScore;scores(i, j, 2)=SSIMScore;
            scores(i, j, 3)=SSIM_L;scores(i, j, 4)=SSIM_S;scores(i, j, 5)=SSIM_R;
        end
    end
end

function p = psnr(x,y)

% psnr - compute the Peack Signal to Noise Ratio, defined by :
%       PSNR(x,y) = 10*log10( max(max(x),max(y))^2 / |x-y|^2 ).
%
%   p = psnr(x,y);
%
%   Copyright (c) 2004 Gabriel Peyr

d = mean( mean( (x(:)-y(:)).^2 ) );
m1 = max( abs(x(:)) );
m2 = max( abs(y(:)) );
m = max(m1,m2);

p = 10*log10( m^2/d );
end

function [mssim, ssim_map, LCS_map] = my_ssim(img1, img2, K, window, L, power)

% ========================================================================
% SSIM Index with automatic downsampling, Version 1.0
% Copyright(c) 2009 Zhou Wang
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereby
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for calculating the
% Structural SIMilarity (SSIM) index between two images
%
% Please refer to the following paper and the website with suggested usage
%
% Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
% quality assessment: From error visibility to structural similarity,"
% IEEE Transactios on Image Processing, vol. 13, no. 4, pp. 600-612,
% Apr. 2004.
%
% http://www.ece.uwaterloo.ca/~z70wang/research/ssim/
%
% Note: This program is different from ssim_index.m, where no automatic
% downsampling is performed. (downsampling was done in the above paper
% and was described as suggested usage in the above website.)
%
% Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (4) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (5) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size
%            depends on the window size and the downsampling factor.
%
%Basic Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim, ssim_map] = ssim(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim, ssim_map] = ssim(img1, img2, K, window, L);
%
%Visualize the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%========================================================================

[M N] = size(img1);
if isempty(K)
    K=[0.01 0.03 0.03];
end
if isempty(window)
    window=fspecial('gaussian', 11, 1.5);
end
if isempty(L)
    L=255;
end
if isempty(power)
    power=[1 1 1];
end

img1 = double(img1);
img2 = double(img2);

% automatic downsampling
f = max(1,round(min(M,N)/256));
%downsampling by f
%use a simple low-pass filter 
if(f>1)
    lpf = ones(f,f);
    lpf = lpf/sum(lpf(:));
    img1 = imfilter(img1,lpf,'symmetric','same');
    img2 = imfilter(img2,lpf,'symmetric','same');

    img1 = img1(1:f:end,1:f:end);
    img2 = img2(1:f:end,1:f:end);
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
C3 = (K(3)*L)^2;
window = window/sum(sum(window));

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 && C2 > 0)   
    L=(2*mu1_mu2+C1)./(mu1_sq+mu2_sq+C1);
    C=(2*sigma12+C2)./(sigma1_sq + sigma2_sq + C2);
    S=(sigma12+C3)./(sqrt(sigma1_sq.*sigma2_sq)+C3);
    LCS_map=zeros(size(L, 1), size(L, 2), 3);
    LCS_map(:, :, 1)=L;LCS_map(:, :, 2)=C;LCS_map(:, :, 3)=S;
    ssim_map=(L.^power(1)).*(C.^power(2)).*(S.^power(3));
 %   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
    fprintf('Not implemented');
    %{
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
	denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
    %}
end

mssim = mean2(ssim_map);
end