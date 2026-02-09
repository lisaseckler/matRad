function [LETmaskDirty,mLET,dij] = matRad_calcLETmask(dij)
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
    if dij.numOfModalities == 2 %exist('dij.original_Dijs','var')
        [i,j,v] = find(dij.original_Dijs{1,k}.physicalDose{1});
        idx = sub2ind(size(dij.original_Dijs{1,k}.physicalDose{1}),i,j);

        if ~isfield(dij.original_Dijs{1,k},'mLETDose')
            %dij.original_Dijs{1,k}.mLETDose{1} = zeros(size(dij.original_Dijs{1,k}.physicalDose{1}));
            dij.original_Dijs{1,k}.LET{1} = (ones(size(dij.original_Dijs{1,k}.physicalDose{1})) * 0.3);
            %dij.original_Dijs{1,k}.mLETDose{1} = dij.original_Dijs{1,k}.LET{1}(idx) .* v;
            %dij.original_Dijs{1,k}.mLETDose{1} = 0.3 * v;
            subIxDirty{1,k} = dij.original_Dijs{1,k}.LET{1} > dij.dirtyDoseThreshold;
            % subIxClean{1,k} = ~subIxDirty{1,k};
        else
            mLET{1,k} = full(dij.original_Dijs{1,k}.mLETDose{1}(idx) ./ v);

            subIxDirty{1,k} = mLET{1,k} > dij.dirtyDoseThreshold;
            % subIxClean{1,k} = ~subIxDirty{1,k};

            mLET{1,k} = sparse(i,j,mLET{1,k},dij.doseGrid.numOfVoxels,dij.original_Dijs{1,k}.totalNumOfBixels);
        end
        LETmaskDirty{1,k} = sparse(i(subIxDirty{1,k}),j(subIxDirty{1,k}),true,dij.doseGrid.numOfVoxels,dij.original_Dijs{1,k}.totalNumOfBixels);
        % LETmaskClean{1,k} = sparse(i(subIxClean{1,k}),j(subIxClean{1,k}),true,dij.doseGrid.numOfVoxels,dij.original_Dijs{1,k}.totalNumOfBixels);
    elseif dij.numOfModalities == 1
        [i,j,v] = find(dij.physicalDose{1});
        idx = sub2ind(size(dij.physicalDose{1}),i,j);

        if ~isfield(dij,'mLETDose')
            dij.LET{1} = (ones(size(dij.physicalDose{1})) * 0.3);
            subIxDirty{1} = dij.LET{1} > dij.dirtyDoseThreshold;
        else
            mLET{1} = full(dij.mLETDose{1}(idx) ./ v);

            subIxDirty{1} = mLET{1} > dij.dirtyDoseThreshold;

            mLET{1} = sparse(i,j,mLET{1},dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
        end
        LETmaskDirty{1} = sparse(i(subIxDirty{1}),j(subIxDirty{1}),true,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
    else
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('No dij found');
    end
end

end