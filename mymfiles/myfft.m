function xx = myfft(x);

[m,n] = size(x);
xx = 1/sqrt(m)*fft(x);
xx = [xx(floor((m+1)/2)+1:end, :); xx(1:floor((m+1)/2), :)];
