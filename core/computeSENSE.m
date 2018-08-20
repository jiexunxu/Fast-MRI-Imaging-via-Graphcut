function x0=computeSENSE(data, sensit, G, lambda, prior)
    conjSensit=conj(sensit);
    b=data;b=Htx(b, conjSensit, G);b=b(:);
    lambda=lambda*lambda;
    y=sensit;
    x0Size=[size(sensit, 1) size(sensit, 2) size(sensit, 3)];
    numPixels=size(sensit, 1)*size(sensit, 2);
    
    if strcmpi(prior, 'TV')
        sizeX=[size(data, 1), size(data, 2), size(data, 3)];
        [nbrList, ~, nbrTermWeights]=computeNbrList(true(sizeX), false(sizeX), [1 1]);
        N=size(nbrList, 1);
        M=sparse([1:N 1:N]', nbrList(:), [nbrTermWeights; -nbrTermWeights], N, sizeX(1)*sizeX(2)*sizeX(3));
        M=sqrt(lambda)*M;Mt=M';MtM=lambda*(M'*M);
        D=numel(data);
      %  x0=pcg(@(x) HtHTVregx(x), b, 1e-6, 80);
        b=zeros(D+size(M, 1), 1);b(1:D)=data(:);
        x0=lsqr(@HxTV, b, [], 80);
    else
        x0=pcg(@(x) HtHregx(x), b, 1e-6, 80);
    end
    x0=reshape(x0, x0Size); 
    % Computes y=(H'H+lambda^2*I)x for multiple coils and multiple slices. Both x and y
    % are in image space
    % Note: y is in the IMAGE SPACE, not k-space
    function z=HtHregx(x)
        y=sensit;
        y=bsxfun(@times, y, reshape(x, x0Size));
        y=fft(fft(y, [], 1), [], 2);
        y=bsxfun(@times, y, G);
        y=conj(fft(fft(conj(y), [], 1), [], 2))/numPixels;
        y=y.*conjSensit;
        z=sum(y, 4);
        z=z(:)+lambda*x;
    end

    % Computes y=(H'H+lambda^2*TV)x for multiple coils and multiple slices. Both x and y
    % are in image space
    % Note: y is in the IMAGE SPACE, not k-space
    function z=HtHTVregx(x)
        y=sensit;
        y=bsxfun(@times, y, reshape(x, x0Size));
        y=fft(fft(y, [], 1), [], 2);
        y=bsxfun(@times, y, G);
        y=conj(fft(fft(conj(y), [], 1), [], 2))/numPixels;
        y=y.*conjSensit;
        z=sum(y, 4);
        z=z(:)+MtM*x;
    end

    function z=HxTV(x, mode)
        if strcmpi(mode, 'transp')
            x1=reshape(x(1:D), size(data));x2=x(D+1:length(x));
            x1=fft(x1, [], 1);
            x1=fft(x1, [], 2);
            x1=bsxfun(@times, x1, G); 
            x1=ifft(x1, [], 1);
            x1=ifft(x1, [], 2);
            x1=x1.*conjSensit;
            x1=sum(x1, 4);
            x2=Mt*x2;
            z=x1(:)+x2;
        else
            z=zeros(D+size(M, 1), 1);
            x=reshape(x, size(data, 1), size(data, 2), size(data, 3));
            y=sensit;
            y=bsxfun(@times, y, x);
            y=fft(y, [], 1);
            y=fft(y, [], 2);
            y=bsxfun(@times, y, G);
            y=ifft(y, [], 1);
            y=ifft(y, [], 2);
            z(1:D)=y(:);
            z(D+1:length(z))=M*x(:);
        end
    end
end