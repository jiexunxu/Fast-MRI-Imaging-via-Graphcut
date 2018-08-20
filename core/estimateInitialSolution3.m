% Compute an initial estimate by finding the analytical solution to
% ||y-Ex||_2^2+lambda*TV(x)_2^2 s.t |Wx|_1<tau
% lambda parameter can be found in the line 
% mu = 0.001;      
function x0=estimateInitialSolution3(data, sensit, G, sigma, lambda, maxIter)
    sizeX=[size(data, 1) size(data, 2) size(data, 3)];   
    [nbrList, ~, nbrTermWeights]=computeNbrList(true(sizeX), false(sizeX), [1 1]);
    N=size(nbrList, 1);
    fullNbrMat=sparse([1:N 1:N]', nbrList(:), [nbrTermWeights; -nbrTermWeights], N, sizeX(1)*sizeX(2)*sizeX(3));
 
    conjSensit=conj(sensit);    
    data=data/max(abs(data(:)));
    wname='db4';decompLvl=6;
    M=numel(data);
    b=zeros(M+size(fullNbrMat, 1), 1);b(1:M)=data(:);
    b=double(b(:));
    fullNbrMat=lambda*fullNbrMat;  
    fullNbrMatTransp=fullNbrMat';
    
    [~, S]=wavedec2(zeros(size(G)), decompLvl, wname);
    opts=spgSetParms('iscomplex', 1, 'iterations', maxIter);    
  %  x0Wavelet=spg_bpdn(@HxTV, b, norm(b)*1e-1, opts);
    x0Wavelet=spgl1(@HxTV, b, 0, norm(b)*sigma, [], opts);
    x0=waverec2(x0Wavelet.', S, wname);    
    
    function y=HxTV(x, mode)
        if mode==1
            % Inverse wavelet transform
            x=x.';
            xreal=waverec2(real(x), S, wname); 
            ximag=waverec2(imag(x), S, wname);
            x=xreal+1i*ximag;
            y=zeros(length(b), 1);
            y(1:M)=Hx(x, G, sensit);
            y(M+1:length(y))=fullNbrMat*x(:);
        else
            xtmp=Htx(reshape(x(1:M), size(data)), conjSensit, G);
            xtmp=xtmp(:)+fullNbrMatTransp*x(M+1:length(x));
            xtmp=reshape(xtmp, size(G));
            [yreal, ~]=wavedec2(real(xtmp), decompLvl, wname);
            [yimag, ~]=wavedec2(imag(xtmp), decompLvl, wname);
            y=yreal+1i*yimag;y=y.';
        end
    end
end