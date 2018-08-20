function xnew = smooth_kspace(x)

method = 1;
szx = size(x);
ii = sqrt(-1);
numrows = szx(1); numcols = szx(2); numslices = szx(3);
hh = [1; 1; 2; 4; 2; 1; 1]; hh = hh/sum(hh);
hpf = [-1; -1; 1; 10; -1; -1; -1]/6;

if method==1
    for sl = 1:numslices
        q = fft2(fftshift(x(:,:,sl)));
        qnew = conv2(abs(q), hh, 'same');
%         for j = 1:numcols
%             qnew(:,j) = cconv(hh, abs(q(:,j)), numrows);
%         end
        qnew = qnew.*exp(ii*angle(q));
        qnew([1,2,end-1,end],:) = q([1,2,end-1,end],:);
        qnew(:,[1,2,end-1,end]) = q(:,[1,2,end-1,end]);      
        xnew(:,:,sl) = fftshift(ifft2(qnew));
    end
end 
