function [E, ens]=smoothnessEnergy(x0, Es_threshold, power, nbrList, varargin)
    if isempty(nbrList)
        E=0;ens=0;
    else
        ens=min(abs(x0(nbrList(:, 1))-x0(nbrList(:, 2))).^power, Es_threshold);
        if nargin==5
            ens=ens.*varargin{1};
        end
        E=sum(ens);    
    end
end