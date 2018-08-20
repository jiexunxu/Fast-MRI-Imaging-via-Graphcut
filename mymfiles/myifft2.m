function xx = myifft2(x);

[m,n] = size(x);
% mx = max(max(abs(x)));
% [I, J] = find(abs(x)==mx);
% if I<m/10 | I>9*m/10 | J<n/10 | J>9*n/10    % check if energy concentrated at corners
	x = [x(:, floor(n/2)+1:end), x(:, 1:floor(n/2))];
	x = [x(floor(m/2)+1:end, :); x(1:floor(m/2), :)];
% end;
xx = sqrt(m*n)*ifft2(x);
