function Hsp = copy_to_sparse(H)
% copies non-sparse ND matrix into sparse 1D matrix
% keeps it memory efficient; meant for non-double non-sparse arrays

szH = size(H);
nr = szH(1);
Hsp = sparse(prod(szH), 1);

for i = 1:nr
    if ndims(H)==2
        q = squeeze(full(H(i,:)));
        Hsp(i,:) = sparse(q);
    elseif ndims(H)==3
        q = squeeze(full(H(i,:,:)));
        Hsp(i,:,:) = sparse(q);
    elseif ndims(H)==4
        q = squeeze(full(H(i,:,:,:)));
        Hsp(i,:,:,:) = sparse(q);
    end
end

       
    