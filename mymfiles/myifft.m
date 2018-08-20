function xx = myifft(x);

[m,n] = size(x);
x = [x(floor(m/2)+1:end, :); x(1:floor(m/2), :)];
xx = sqrt(m)*ifft(x);
