function [ Data_Dose,Data_Fluence,Zdepths, rDepths, AreaOfScoringRings] = combineScoringVolumes(Data_In, Data_Fluence_In,resolution, R0)

    R0 = resolution.zCoarseNum + 0.25*R0;
    R = 0:resolution.r:resolution.r*resolution.rNum;
    AreaOfScoringRings = pi.*(R(2:end).^2-R(1:end-1).^2);
    % Middle of scoring Voxel
    ZdepthsCoarse = resolution.zCoarse/2 : resolution.zCoarse : (resolution.zCoarseNum-0.5)*resolution.zCoarse;
    ZdepthsFine = R0-10+resolution.zFine/2 : resolution.zFine : R0+10-resolution.zFine/2;
    rDepths = resolution.r/2:resolution.r : (resolution.rNum-0.5)*resolution.r;

    % Combine coarse and fine zScoring
    IxStart = find((ZdepthsFine(1) - ZdepthsCoarse)<=0);
    IxStart = IxStart(1);
    IxEnd = find((ZdepthsFine(end) - ZdepthsCoarse)>=0);
    IxEnd = IxEnd(end);

    Zdepths = ZdepthsCoarse(1:IxStart-1);
    Zdepths = [Zdepths, ZdepthsFine];
    Zdepths = [Zdepths, ZdepthsCoarse(IxEnd+1:end)];

    Data_Dose.ED = combine(Data_In.ED_Coarse,Data_In.ED_Fine,IxStart,IxEnd);

    Data_Fluence.energyBin = Data_Fluence_In.energyBin;
    Data_Fluence.energyBinEl = Data_Fluence_In.energyBinEl; 
    
    for i=1:length(Data_Fluence_In.Z)
        Data_Fluence.spectra(i).Z = Data_Fluence_In.Z(i);
        Data_Fluence.spectra(i).A = Data_Fluence_In.A(i);
        Data_Fluence.spectra(i).fluenceSpectrum = combine(Data_Fluence_In.spectrum(i).Coarse', Data_Fluence_In.spectrum(i).Fine', IxStart,IxEnd)'*(pi*R(end)*R(end));
        
    end


end

function result = combine(InCoarse,InFine,IxStart,IxEnd)
   result  = InCoarse(1:IxStart-1,:);
   result = [result; InFine];
   result = [result; InCoarse(IxEnd+1:end,:)];
end