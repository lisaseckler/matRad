function dij = matRad_calcDirtyDose(LET_thres,dij)
% Calculates Dirty and Clean Dose by using LET threshold
%
% call
%   dij = matRad_calcDirtyDose(LET_thres,dij)
%
% input
%   LET_thres:  LET threshold, above: dirty dose, below: clean dose
%   dij:        matRad dij struct
%
% output
%   dij:        matRad dij struct with dirty and clean dose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('dij','var') || isempty(dij)
    disp('Not enough input arguments! Calculation is not working.')
else
    if ~exist('LET_thres','var') || isempty(LET_thres)
        disp('Not enough input arguments! Calculation is not working.')
    else
        dij.dirtyDoseThreshold                       = LET_thres;
        [dij.LETmaskDirty,dij.LETmaskClean,dij.mLET] = matRad_calcLETmask(dij);

        for k = 1:dij.numOfModalities
            dij.dirtyDose{1,k}                       = dij.LETmaskDirty{1,k} .* dij.original_Dijs{1,k}.physicalDose{1};
            dij.cleanDose{1,k}                       = dij.LETmaskClean{1,k} .* dij.original_Dijs{1,k}.physicalDose{1};
        end
        
        DirtySize = zeros(size(dij.original_Dijs{1,1}.physicalDose{1}));
        LETDirty_Size = size(dij.dirtyDose{1,2});
        row = LETDirty_Size(1,1);
        column = LETDirty_Size(1,2);
        DirtySize(1:row,1:column) = dij.dirtyDose{1,2};
        photonDD = sparse(DirtySize);
        dij.dirtyDose{1,3} = dij.dirtyDose{1,1} + photonDD;
        DirtySize(1:row,1:column) = dij.cleanDose{1,2};
        photonCD = sparse(DirtySize);
        dij.cleanDose{1,3} = dij.cleanDose{1,1} + photonCD;

    end
end

end