function B=waveletTransform3D(A, w1, w2, direction)
    alpha=sum(w1);w1=w1/alpha;w2=w2/alpha;
    if direction==+1
        B=waveletDecomposition3D(A, w1, w2);
    elseif direction==-1
        B=waveletReconstruction3D(A, w1, w2);
    else
        B=0;
    end
end

% LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH
function B=waveletDecomposition3D(A, Lo_D, Hi_D)
    [m, n, p]=size(A);
    q=length(Lo_D);
    B=complex(zeros(m+q-1, n+q-1, (p+q-1)*8), 0);    
    BL=directionalConv(A, Lo_D, 1);
    BH=directionalConv(A, Hi_D, 1);
    BLL=directionalConv(BL, Lo_D, 2);
    BLH=directionalConv(BL, Hi_D, 2);
    BHL=directionalConv(BH, Lo_D, 2);
    BHH=directionalConv(BH, Hi_D, 2);
    L=p+q-1;
    B(:, :, 1:L)=directionalConv(BLL, Lo_D, 3);
    B(:, :, L+1:2*L)=directionalConv(BLL, Hi_D, 3);
    B(:, :, 2*L+1:3*L)=directionalConv(BLH, Lo_D, 3);
    B(:, :, 3*L+1:4*L)=directionalConv(BLH, Hi_D, 3);
    B(:, :, 4*L+1:5*L)=directionalConv(BHL, Lo_D, 3);
    B(:, :, 5*L+1:6*L)=directionalConv(BHL, Hi_D, 3);
    B(:, :, 6*L+1:7*L)=directionalConv(BHH, Lo_D, 3);
    B(:, :, 7*L+1:8*L)=directionalConv(BHH, Hi_D, 3);
end

function B=waveletReconstruction3D(A, Lo_R, Hi_R)
    [m, n, p]=size(A);p=p/8;
    q=length(Lo_R);
    L=p+q-1;
    B=zeros(m+q-1, n+q-1, p+q-1);    
    B=B+singleBlockRecon(A(:, :, 1:L), Lo_R, Lo_R, Lo_R);
    B=B+singleBlockRecon(A(:, :, L+1:2*L), Hi_R, Lo_R, Lo_R);
    B=B+singleBlockRecon(A(:, :, 2*L+1:3*L), Lo_R, Hi_R, Lo_R);
    B=B+singleBlockRecon(A(:, :, 3*L+1:4*L), Hi_R, Hi_R, Lo_R);
    B=B+singleBlockRecon(A(:, :, 4*L+1:5*L), Lo_R, Lo_R, Hi_R);
    B=B+singleBlockRecon(A(:, :, 5*L+1:6*L), Hi_R, Lo_R, Hi_R);
    B=B+singleBlockRecon(A(:, :, 6*L+1:7*L), Lo_R, Hi_R, Hi_R);
    B=B+singleBlockRecon(A(:, :, 7*L+1:8*L), Hi_R, Hi_R, Hi_R);    
    B=B(q:m, q:n, q:p);
end

function B=directionalConv(A, w, dir)
    [m, n, p]=size(A);
    q=length(w);
    if dir==1   
        B=complex(zeros(m+q-1, n, p), 0);
        for j=1:n
            for k=1:p
                v=conv(squeeze(A(:, j, k)), w);
                B(:, j, k)=v;
            end
        end
    elseif dir==2
        B=complex(zeros(m, n+q-1, p), 0);
        for i=1:m
            for k=1:p
                v=conv(squeeze(A(i, :, k)), w);
                B(i, :, k)=v;
            end
        end
    elseif dir==3 
        B=complex(zeros(m, n, p+q-1), 0);
        for i=1:m
            for j=1:n
                v=conv(squeeze(A(i, j, :)), w);
                B(i, j, :)=v;
            end
        end
    end
end

% Conv in z, y, x direction with w1, w2, w3
function B=singleBlockRecon(A, w1, w2, w3)
    [m, n, p]=size(A);
    q=length(w1);
    B=complex(zeros(m+q-1, n+q-1, p+q-1), 0);
    for i=1:m
        for j=1:n
            v=conv(squeeze(A(i, j, :)), w1);
            B(i, j, :)=v;
        end
    end
    
    for i=1:m
        for k=1:p+q-1
            v=conv(squeeze(B(i, 1:n, k)), w2);
            B(i, :, k)=v;
        end
    end
    
    for j=1:n+q-1
        for k=1:p+q-1
            v=conv(squeeze(B(1:m, j, k)), w3);
            B(:, j, k)=v;
        end
    end
end