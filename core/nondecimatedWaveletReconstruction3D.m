function B=nondecimatedWaveletReconstruction3D(A)
    LR=[2.30377813308855e-001 7.14846570552542e-001 6.30880767929590e-001 -2.79837694169838e-002 -1.87034811718881e-001 3.08413818359870e-002 3.28830116669829e-002 -1.05974017849973e-002];
    HR=[-1.05974017849973e-002 -3.28830116669829e-002 3.08413818359870e-002 1.87034811718881e-001 -2.79837694169838e-002 -6.30880767929590e-001 7.14846570552542e-001 -2.30377813308855e-001];
    alpha=sum(LR);LR=LR/alpha;HR=HR/alpha;
    B=waveletReconstruction3D(A, LR, HR);
end

function B=waveletReconstruction3D(A, Lo_R, Hi_R)
    [m, n, p]=size(A);p=p/8;
    q=length(Lo_R);L=size(A, 3)/8;
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