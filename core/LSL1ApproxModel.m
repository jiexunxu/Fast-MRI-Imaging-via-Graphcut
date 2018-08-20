% Approximates objective with ||H(x0+V*alpha)-y||_2^2+EsEdratio*||DM(x0+V*alpha)||_1
function coefs=LSL1ApproxModel(V, x0, sensitivities, G, data, objParams, nbrTermWeights, fullNbrMat, EsEdratio, LSIter, smoothScaling)
    b1=data(:)-Hx(x0, G, sensitivities);Htb1=Htx(b1, sensitivities, G);
    DMV=spdiags(nbrTermWeights, 0, length(nbrTermWeights), length(nbrTermWeights))*fullNbrMat*V;
    b2=spdiags(nbrTermWeights, 0, length(nbrTermWeights), length(nbrTermWeights))*fullNbrMat*x0(:);
    Ed=norm(Hx(reshape(V*x0(:), size(x0)), G, sensitivities)-b1, 2)^2;
    Es=norm(DMV*x0(:)+b2, 1);
    EsEdratio=double(EsEdratio*Ed/Es*smoothScaling);
    opts=optimset('GradObj', @LSL1ApproxFunGrad, 'MaxIter', 30);
    coefs=fminunc(@LSL1ApproxFun, x0(:), opts);
    
    
    function y=LSL1ApproxFun(x)
        y=Hx(reshape(V*x, size(x0)), G, sensitivities)-b1;
        y=sum(y.*y);
        y=y+EsEdratio*norm(DMV*x+b2, 1);
    end

    function grad=LSL1ApproxFunGrad(x)
        grad=HtHx(reshape(V*x, size(x0)), sensitivities, G)-2*Htb1;
        signs=sign(DMV*x+b2);
        grad=grad+EsEdratio*(DMV'*signs);
    end
end