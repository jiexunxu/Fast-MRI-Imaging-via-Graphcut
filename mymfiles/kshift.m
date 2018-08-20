function AA = kshift(A, ind1, ind2)
% centers A in k-space so that peak k value is in the centre
% cyclically shifts k-space - prob need to change?
% you can supply offset values in row and col coords (ind1, ind2),
% or let kshift do it for you (dont pass ind1 or ind2 if so)
% supply only one offset if row_offset = col_offset

q = A;
% lpf = [0 1 0; 1 2 1; 0 1 0];  lpf = lpf/sum(sum(lpf));
% q = conv2(A, lpf, 'same');
[m,n] = size(q);

if nargin == 3  dr = ind1;  dc = ind2;
elseif nargin==2  dr = ind1(1);  dc = ind1(2);
else  
    mx = max(max(abs(q)));
    [centrer, centrec] = find(abs(q)==mx);
    dr = round(centrer-m/2);
    dc = round(centrec-n/2);
end
q = A;
    
if dr<0
	q  = [q(m-abs(dr)+1:m, :); q(1:m-abs(dr), :)];
elseif dr>0
	q  = [q(abs(dr)+1:m, :); q(1:abs(dr), :)];
end
[m,n] = size(q);
if dc<0
    q = [q(:, end-abs(dc)+1:end), q(:, 1:end-abs(dc))];
elseif dc>0
    q = [q(:, abs(dc)+1:end), q(:, 1:abs(dc))];
end
AA = q;
