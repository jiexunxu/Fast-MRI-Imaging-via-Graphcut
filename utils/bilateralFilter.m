function x0blurred=bilateralFilter(x0, w, sigma)
    x0blurred=zeros(size(x0));
    for i=1:size(x0, 3) 
        img=abs(x0(:, :, i));maxClr=max(abs(img(:)));             
        x0blurred(:, :, i)=bfilter2(img/maxClr, w, sigma);
        x0blurred(:, :, i)=x0blurred(:, :, i)*maxClr;
    end
    x0blurred=x0blurred(:).*exp(1i*angle(x0(:)));
    x0blurred=reshape(x0blurred, size(x0));
end