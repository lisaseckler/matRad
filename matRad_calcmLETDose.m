function dij = matRad_calcmLETDose(dij,pln)
% calculating mLETDose using mixed modalities
%
% call
%   dij = matRad_calcmLETDose(dij,pln)
%
% input
%   pln:    matRad pln struct
%   dij:    matRad dij struct
%
% output
%   dij:        matRad dij struct with dirty and clean dose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dij.mLETDose = {};
 m = 2;
        rad = ["protons","photons","carbon","helium"];
        for k = 1:dij.numOfModalities
            if ~isfield(dij.original_Dijs{1,k},'mLETDose')
                dij.original_Dijs{1,k}.mLETDose{1} = zeros(size(dij.original_Dijs{1,k}.physicalDose{1}));
            end
        end

        if pln(2).radiationMode == rad(1,2)
            Size = zeros(size(dij.original_Dijs{1,m-1}.physicalDose{1}));
            LET_Size = size(dij.original_Dijs{1,m}.mLETDose{1});
            row = LET_Size(1,m-1);
            column = LET_Size(1,m);
            Size(1:row,1:column) = dij.original_Dijs{1,m}.mLETDose{1};
            photonmLETDose = sparse(Size);
            dij.mLETDose{1,1} = dij.original_Dijs{1,m-1}.mLETDose{1};
            dij.mLETDose{1,2} = dij.original_Dijs{1,m}.mLETDose{1};
            dij.mLETDose{1,3} = dij.mLETDose{1,1} + photonmLETDose;
            
        elseif pln(2).radiationMode == rad(1,3)
            Size = zeros(size(dij.original_Dijs{1,m}.physicalDose{1}));
            LET_Size = size(dij.original_Dijs{1,m-1}.mLETDose{1});
            row = LET_Size(1,1);
            column = LET_Size(1,2);
            Size(1:row,1:column) = dij.original_Dijs{1,m}.mLETDose{1};
            protonmLETDose = sparse(Size);
            dij.mLETDose{1,1} = dij.original_Dijs{1,m-1}.mLETDose{1};
            dij.mLETDose{1,2} = dij.original_Dijs{1,m}.mLETDose{1};
            dij.mLETDose{1,3} = dij.mLETDose{1,2} + protonmLETDose;
        end

end