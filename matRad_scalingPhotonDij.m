function dij = matRad_scalingPhotonDij(dij)
% preconditioning: scaling photon dij by 100
dij.original_Dijs{1,2}.physicalDose{1} = dij.original_Dijs{1,2}.physicalDose{1}./100;
dij.original_Dijs{1,2}.dirtyDose{1} = dij.original_Dijs{1,2}.dirtyDose{1}./100;
dij.original_Dijs{1,2}.mLETDose{1} = dij.original_Dijs{1,2}.mLETDose{1}./100;
end