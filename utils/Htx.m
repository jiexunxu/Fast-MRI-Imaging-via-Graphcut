% Computes y=H'*y
function y=Htx(y, conjSensit, G)
    y=fft(y, [], 1);
    y=fft(y, [], 2);
    y=bsxfun(@times, y, G);
 %   y=fft(conj(y), [], 1);
 %   y=fft(y, [], 2);
 %   y=conj(y)/(size(y, 1)*size(y, 2));  
    y=ifft(y, [], 1);
    y=ifft(y, [], 2);

    y=y.*conjSensit;
    y=sum(y, 4);
end