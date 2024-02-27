function dij = matRad_calcDirtyDose(LET_thres,dij,pln)
% Calculates Dirty and Clean Dose by using LET threshold
%
% call
%   dij = matRad_calcDirtyDose(LET_thres,dij)
%
% input
%   LET_thres:  LET threshold, above: dirty dose, below: clean dosedij
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
        [dij.LETmaskDirty,dij.mLET,~] = matRad_calcLETmask(dij);

        for k = 1:dij.numOfModalities
            dij.original_Dijs{1,k}.dirtyDoseThreshold    = LET_thres;
            dij.original_Dijs{1,k}.dirtyDose{1}             = dij.LETmaskDirty{1,k} .* dij.original_Dijs{1,k}.physicalDose{1};
            % dij.original_Dijs{1,k}.cleanDose             = dij.LETmaskClean{1,k} .* dij.original_Dijs{1,k}.physicalDose{1};
        end
        
        % m = 2;
        % rad = ["protons","photons","carbon","helium"];

        % for adding up the dirty dose
        % if dij.numOfModalities == 2
        %     if pln(2).radiationMode == rad(1,2)
        %         DirtySize = zeros(size(dij.original_Dijs{1,m-1}.physicalDose{1}));
        %         LETDirty_Size = size(dij.original_Dijs{1,2}.dirtyDose{1});
        %         row = LETDirty_Size(1,m-1);
        %         column = LETDirty_Size(1,m);
        %         DirtySize(1:row,1:column) = dij.original_Dijs{1,2}.dirtyDose{1};
        %         photonDD = sparse(DirtySize);
        %         dij.dirtyDose{1} = dij.original_Dijs{1,1}.dirtyDose{1} + photonDD;
        %         % DirtySize(1:row,1:column) = dij.original_Dijs{1,2}.cleanDose;
        %         % photonCD = sparse(DirtySize);
        %         % dij.cleanDose = dij.original_Dijs{1,1}.cleanDose + photonCD;
        %     elseif pln(2).radiationMode == rad(1,3)
        %         DirtySize = zeros(size(dij.original_Dijs{1,m}.physicalDose{1}));
        %         LETDirty_Size = size(dij.original_Dijs{1,1}.dirtyDose{1});
        %         row = LETDirty_Size(1,m-1);
        %         column = LETDirty_Size(1,m);
        %         DirtySize(1:row,1:column) = dij.original_Dijs{1,1}.dirtyDose{1};
        %         protonDD = sparse(DirtySize);
        %         dij.dirtyDose{1} = dij.original_Dijs{1,2}.dirtyDose + protonDD;
        %         % DirtySize(1:row,1:column) = dij.original_Dijs{1,1}.cleanDose;
        %         % protonCD = sparse(DirtySize);
        %         % dij.cleanDose = dij.original_Dijs{1,2}.cleanDose + protonCD;
        %     end
        % end

    end
end

end