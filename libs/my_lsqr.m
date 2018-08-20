function [x,flag,relres,iter,resvec,lsvec] = lsqr(afun,b,n,tol,maxit,varargin)
if isempty(tol)
    tol=1e-6;
end
if isempty(maxit)
    maxit=20;
end
% Set up for the method
n2b = norm(b);                     % Norm of rhs vector, b
flag = 1;
tolb = tol * n2b;                  % Relative tolerance
u = b;

beta = norm(u);
% Norm of residual r=b-A*x is estimated well by prod_i abs(sin_i)
normr = beta;
if beta ~= 0
    u = u / beta;
end
c = 1;
s = 0;
phibar = beta;
v = iterapp('mtimes',afun,atype,afcnstr,u,varargin{:},'transp');
alpha = norm(v);
if alpha ~= 0
    v = v / alpha;
end
d = zeros(n,1);

% norm((A*inv(M))'*r) = alpha_i * abs(sin_i * phi_i)
normar = alpha * beta;

x = zeros(n,1);
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
stag = 0;                      % stagnation of the method
iter = maxit;                  % Assume lack of convergence until it happens
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit  
    u = iterapp('mtimes',afun,atype,afcnstr,v,varargin{:},'notransp') - alpha * u;
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
    if (phi == 0)              % stagnation of the method
        stag = 1;
    end
    phibar = s * phibar;
    d = (z - thet * d) / rho;
    sumnormd2 = sumnormd2 + (norm(d))^2;
    
    % Check for stagnation of the method
    if abs(phi)*norm(d) < eps*norm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    if stag >= maxstagsteps
        flag = 3;
        iter = ii-1;
        resvec = resvec(1:iter+1);
        lsvec = lsvec(1:iter+1);
        break
    end
    
    x = x + phi * d;
    normr = abs(s) * normr;
    resvec(ii+1) = normr;
    vt = iterapp('mtimes',afun,atype,afcnstr,u,varargin{:},'transp');    
    v = vt - beta * v;
    alpha = norm(v);
    v = v / alpha;
    normar = alpha * abs( s * phi);
    
end                            % for ii = 1 : maxit

if flag == 1
    if normar/(norma*normr) <= tol % check for convergence in min{|b-A*x|}
        flag = 0;
        iter = maxit;
    end
    
    if normr <= tolb           % check for convergence in A*x=b
        flag = 0;
        iter = maxit;
    end
end

relres = normr/n2b;

% only display a message if the output flag is not used
if nargout < 2
    itermsg('lsqr',tol,maxit,ii,flag,iter,relres);
end
