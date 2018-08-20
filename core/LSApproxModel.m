% Approximates objective with ||H(x0+V*alpha)-y||_2^2+EsEdratio*||DM(x0+V*alpha)||_2^2
function coefs=LSApproxModel(V, x0, sensitivities, conjSensitivities, G, data, objParams, nbrTermWeights, fullNbrMat, EsEdratio, LSIter, smoothScaling)
    DMV=spdiags(nbrTermWeights, 0, length(nbrTermWeights), length(nbrTermWeights))*fullNbrMat*V;
    Vt=V';DMVt=DMV';
    z=sensitivities;N=numel(sensitivities);M=N+size(fullNbrMat, 1);
    w1=sqrt(objParams(1));w2=sqrt(objParams(2));
    hx=sensitivities;hx=Hx(x0, G, hx);hx=data(:)-hx(:);
    
    b=[w1*hx; nbrTermWeights.*(fullNbrMat*x0(:))*(-w2)];
    Ed=norm(b(1:numel(sensitivities)));Es=norm(b(numel(sensitivities)+1:length(b)));
    SmoothScale=double(sqrt(EsEdratio)*Ed/Es*smoothScaling);w2=w2*double(SmoothScale);
    b=[w1*hx; nbrTermWeights.*(fullNbrMat*x0(:))*(-w2)];
    b=double(b);DMV=DMV*w2;DMVt=DMVt*w2;
    coefs=my_lsqr(@ApproxLSFun, b, LSIter, zeros(size(V, 2), 1));
 %   coefs=my_lsqr(@ApproxLSFun, b, LSIter, zeros(size(V, 2), 1));
    % || w1*(HV*alpha + (Hx0 - y)) ||+|| w2*(DMV*alpha + DM*x0) || =>
    % || [w1*HV; w2*DMV] * alpha - [w1*(y-Hx0); -w2*DM*x0] ||    
    function y=ApproxLSFun(x, arg)
        if strcmp(arg, 'transp')            
            z=reshape(x(1:N), size(sensitivities));
            z=Htx(z, conjSensitivities, G);
            y=Vt*double(z(:))*w1+DMVt*x(N+1:M);
        else
            z=sensitivities;
            z=Hx(reshape(V*x, size(x0)), G, z);
            y=[w1*double(z(:)); DMV*x];  
        end
        y=double(y);
    end         
end

function x = my_lsqr(afun,b,maxit,x0)
    n=size(x0, 1);    
    % Assign default values to unspecified parameters

    xInit = ~isempty(x0);
    x=x0;

    u = b;
    if xInit
        u = u - afun(x, 'notransp');
    end
    beta = norm(u);
    % Norm of residual r=b-A*x is estimated well by prod_i abs(sin_i)
    normr = beta;
    if beta ~= 0
        u = u / beta;
    end
    c = 1;
    s = 0;
    phibar = beta;
    v = afun(u, 'transp');
    alpha = norm(v);
    if alpha ~= 0
        v = v / alpha;
    end
    d = zeros(n,1);

    % norm((A*inv(M))'*r) = alpha_i * abs(sin_i * phi_i)
    normar = alpha * beta;
    % And make sure we have x if it was not initialized.
    if ~xInit
        x = zeros(n,1);
    end

    % Poorly estimate norm(A*inv(M),'fro') by norm(B_{ii+1,ii},'fro')
    % which is in turn estimated very well by
    % sqrt(sum_i (alpha_i^2 + beta_{ii+1}^2))
    norma = 0;
    % norm(inv(A*inv(M)),'fro') = norm(D,'fro')
    % which is poorly estimated by sqrt(sum_i norm(d_i)^2)
    sumnormd2 = 0;
    resvec = zeros(maxit+1,1);     % Preallocate vector for norm of residuals
    resvec(1) = normr;             % resvec(1,1) = norm(b-A*x0)
    lsvec = zeros(maxit,1);        % Preallocate vector for least squares estimates

% loop over maxit iterations (unless convergence or failure)

    for ii = 1 : maxit
        z=v;
        u = afun(z, 'notransp') - alpha * u;
        beta = norm(u);
        u = u / beta;
        norma = norm([norma alpha beta]);
        lsvec(ii) = normar / norma;
        thet = - s * alpha;
        rhot = c * alpha;
        rho = sqrt(rhot^2 + beta^2);
        c = rhot / rho;
        s = - beta / rho;
        phi = c * phibar;
        phibar = s * phibar;
        d = (z - thet * d) / rho;
        sumnormd2 = sumnormd2 + (norm(d))^2;

        x = x + phi * d;
        normr = abs(s) * normr;
        resvec(ii+1) = normr;
        vt = afun(u, 'transp');
        v = vt - beta * v;
        alpha = norm(v);
        v = v / alpha;
        normar = alpha * abs( s * phi);
    end
end