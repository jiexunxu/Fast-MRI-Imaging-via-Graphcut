function xx = myfft2(x);

[m,n] = size(x);
xx = 1/sqrt(m*n)*fft2(x);
% mx = max(max(abs(xx)));
% [I, J] = find(abs(xx)==mx);
% if I<m/10 | I>9*m/10 | J<n/10 | J>9*n/10    % check if energy concentrated at corners
	xx = [xx(floor((m+1)/2)+1:end, :); xx(1:floor((m+1)/2), :)];
	xx = [xx(:, floor((n+1)/2)+1:end), xx(:, 1:floor((n+1)/2))];
% end;