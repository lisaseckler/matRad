function [LETmaskDirty,LETmaskClean,mLET] = matRad_calcLETmask(dij)
% Calculates logical matrix where LET is above (LETmask.Dirty) or below (LETmask.Clean) a certain threshold
%
% call
%   LETmask = matRad_calcLETmask(dij)
%
% input
%   dij:       matRad dij struct
%
% output
%   LETmask:   logical matrix for dirty dose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% instead of:
% mLET = full(dij.mLETDose{1}(:,:))./full(dij.physicalDose{1}(:,:));
% LETmask = mLET > dij.dirtyDoseThreshold;

% logical matrix is calculated like this:

mLET = {};
subIxDirty = {};
subIxClean = {};
LETmaskDirty = {};
LETmaskClean = {};

for k = 1:dij.numOfModalities
    
    if ~isfield(dij.original_Dijs{1,k},'mLETDose')
        %dij.original_Dijs{1,k}.mLETDose{1} = zeros(size(dij.original_Dijs{1,k}.physicalDose{1}));
        dij.original_Dijs{1,k}.mLETDose{1} = ones(size(dij.original_Dijs{1,k}.physicalDose{1})) * 0.3;
    end

    [i,j,v] = find(dij.original_Dijs{1,k}.physicalDose{1});
    idx = sub2ind(size(dij.original_Dijs{1,k}.physicalDose{1}),i,j);

    mLET{1,k} = full(dij.original_Dijs{1,k}.mLETDose{1}(idx) ./ v);  
    
    subIxDirty{1,k} = mLET{1,k} > dij.dirtyDoseThreshold;
    subIxClean{1,k} = ~subIxDirty{1,k};

    mLET{1,k} = sparse(i,j,mLET{1,k},dij.doseGrid.numOfVoxels,dij.original_Dijs{1,k}.totalNumOfBixels);
    LETmaskDirty{1,k} = sparse(i(subIxDirty{1,k}),j(subIxDirty{1,k}),true,dij.doseGrid.numOfVoxels,dij.original_Dijs{1,k}.totalNumOfBixels);
    LETmaskClean{1,k} = sparse(i(subIxClean{1,k}),j(subIxClean{1,k}),true,dij.doseGrid.numOfVoxels,dij.original_Dijs{1,k}.totalNumOfBixels);
end
    
end