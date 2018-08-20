function x2=contrastEnhance(x, alpha)
    x2=abs(x)/max(abs(x(:)));
    x2=imadjust(x2, alpha, [0 1]);
    x2=x2*max(abs(x(:)));    
    x2=x2.*exp(1i*angle(x));
end