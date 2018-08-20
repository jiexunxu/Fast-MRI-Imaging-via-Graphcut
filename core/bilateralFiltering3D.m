function B=bilateralFiltering3D(A, w, sigD, sigR)    
    Br=bilateralFiltering3D_mex(real(A), size(A, 1), size(A, 2), size(A, 3), w, sigD, sigR);
    Bi=bilateralFiltering3D_mex(imag(A), size(A, 1), size(A, 2), size(A, 3), w, sigD, sigR);
    B=Br+1i*Bi;B=reshape(B, size(A));    
end
%{
function B=bilateralFiltering3D(A, w, sigD, sigR)

    Br=filtering3D(real(A), w, sigD, sigR);
    Bi=filtering3D(imag(A), w, sigD, sigR);
    B=Br+1i*Bi;

end

function B=filtering3D(A, w, sigD, sigR)
    [m, n, p]=size(A);
    B=zeros(m, n, p);
    wx=w(1);wy=w(2);wz=w(3);
    % Gaussian distance weights
    [X, Y, Z]=meshgrid(-wx:wx, -wy:wy, -wz:wz);
    %{
    X=zeros(wx*2+1, wy*2+1, wz*2+1);
    Y=zeros(wx*2+1, wy*2+1, wz*2+1);
    Z=zeros(wx*2+1, wy*2+1, wz*2+1);
    for i=1:wx*2+1
        Y(i, :, :)=i-wx-1;
    end
    for i=1:wy*2+1
        X(:, i, :)=i-wy-1;
    end
    for i=1:wz*2+1
        Z(:, :, i)=i-wz-1;
    end
    %}
    G =exp(-(X.^2+Y.^2+Z.^2)/(2*sigD^2));
    for x=wx+1:m-wx
        for y=wy+1:n-wy
            for z=wz+1:p-wz
                 % Local region
                 I = A(x-wx:x+wx, y-wy:y+wy, z-wz:z+wz);
                 % Intensity weights
                 H=exp(-(I-A(x, y, z)).^2/(2*sigR^2));     
                 % Calculate bilateral filter response.
                 F = H.*G;
                 B(x, y, z) = sum(F(:).*I(:))/sum(F(:));
            end
        end
    end
end
%}