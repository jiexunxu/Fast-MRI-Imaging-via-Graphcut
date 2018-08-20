function B=nondecimatedWaveletTransform3D(A)
    LD=[-1.05974017849973e-002 3.28830116669829e-002 3.08413818359870e-002 -1.87034811718881e-001 -2.79837694169838e-002 6.30880767929590e-001 7.14846570552542e-001 2.30377813308855e-001];
    HD=[-2.30377813308855e-001 7.14846570552542e-001 -6.30880767929590e-001 -2.79837694169838e-002 1.87034811718881e-001 3.08413818359870e-002 -3.28830116669829e-002 -1.05974017849973e-002];
    alpha=sum(LD);LD=LD/alpha;HD=HD/alpha;
    B=waveletDecomposition3D(A, LD, HD);
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