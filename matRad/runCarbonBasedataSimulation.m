%% Introduction
% script that runs my Helium Base Data Simulations
% It creates the corresponding Topas infput file from the Basic input File,
% and adds the current energy, focus index and so on
% Then topas is called and the results are saved

%% Load Energys to Simulate, inital Focus of these Energys and Energy Spread
% load BeamOptics.mat
% load EnergySpectrumData.mat
% load carbon_HITfixedBL.mat
%% for Helium
% alpha = 2.528e-3;
% p = 1.746;
% approxRange = alpha*Eall.^p*10;
% 
% simulateHalfLength = round(0.75*approxRange,0);
approxRange = [];

energyIx = 1:255;

for i = energyIx
    approxRange(i,1) = machine.data(i).peakPos;
end

simulateHalfLength = round(approxRange,0);
%% Simulation Loop

folderName = 'Scripts_big';
workingDir = pwd;
nPrimary = 1000000;
Emax = zeros(1,numel(E));
Ebins = zeros(1,numel(E));



for i = 1:numel(E)

    currentEnergy       = E(i);
    currParams.Ebins    = EnergyBins{i};
    currParams.Espec    = EnergySpectrum{i}/sum(EnergySpectrum{i}); %normaliue
    currParams.HL       = simulateHalfLength(i);
    currParams.range    = approxRange(i);
    currParams.sigmaPrime   = sigmaT(i);
    currParams.correlation      = rho(i);
    currParams.sigma    = sigmaSqIso(i);
    currParams.Emax = round(max(currParams.Ebins) + 30);
    currParams.Ebins_scoring = round(currParams.Emax/12);
    Emax(i) = currParams.Emax;
    Ebins(i) = currParams.Ebins_scoring;
   
    fprintf('Calling TOPAS simulation of Energy %g MeV\n', currentEnergy);
    mkdir([folderName, filesep,  sprintf('Energy%g',currentEnergy)]);
    mkdir([folderName, filesep,  sprintf('Energy%g',currentEnergy), filesep, 'Results']);
    mkdir([folderName, filesep,  sprintf('Energy%g',currentEnergy), filesep, 'Results' filesep 'Fluence']);
    % write TOPAS data base file 
    writeSetupFile(workingDir,folderName, currentEnergy,currParams,nPrimary);
        
    writeCallFile(workingDir,folderName, currentEnergy); 
  
end
   
disp('All Simulations Files Written')
fclose all;

%%
save('SimuParams.mat','nPrimary','simulateHalfLength','approxRange','E','EnergyBins','EnergySpectrum','sigmaT','rho','sigmaSqIso','Emax','Ebins')
%% Function that write the TOPAS simulation file
function writeSetupFile(workingDir, folderName, energy, currParams,nPrimary)

    fileName = [filesep, sprintf('Energy%g',energy), filesep, 'carbon_DoseAndRBE.txt']; 

    fileID = fopen(fullfile(workingDir,folderName, fileName),'w+');
    basicFile = fileread([folderName filesep 'ion_DoseAndRbe_basic.txt']);
    fprintf(fileID,'%s\n\n',basicFile);
    fprintf(fileID, ' \n');

    fprintf(fileID, 'd:Sim/HL = %f mm \n', currParams.HL );
    fprintf(fileID, 'd:Sim/Range = %f mm \n', currParams.range);
    fprintf(fileID, 'i:Sim/ZCoarse/ZBins = %i \n', currParams.HL );
    fprintf(fileID, ' \n');
    
    fprintf(fileID,'s:So/MySource/BeamEnergySpectrumType = "Continuous"\n')
    fprintf(fileID, 'dv:So/MySource/BeamEnergySpectrumValues = %d ', length(currParams.Ebins));
    fprintf(fileID, '%.6f ', currParams.Ebins);
    fprintf(fileID, 'MeV\n');
    fprintf(fileID, 'uv:So/MySource/BeamEnergySpectrumWeights= %d ', length(currParams.Espec));
    fprintf(fileID, '%.6f ', currParams.Espec);
    fprintf(fileID, '\n');
    fprintf(fileID, 'd:Sim/MySource/SigmaX = %f mm \n', currParams.sigma);
    fprintf(fileID, 'u:Sim/MySource/SigmaXprime = %f \n', currParams.sigmaPrime);
    fprintf(fileID, 'u:Sim/MySource/CorrelationX = %f \n', currParams.correlation);
    fprintf(fileID, 'i:Sim/MySource/NumberOfHistories = %i \n', nPrimary)
    fprintf(fileID, '\n');
    fprintf(fileID, 'd:Sim/EMax_ion = %f MeV \n',currParams.Emax );
    fprintf(fileID, 'd:Sim/EMin_ion = %f MeV \n', 0);
    fprintf(fileID, 'i:Sim/Ebins_ion = %i \n',currParams.Ebins_scoring);
    fclose(fileID);

    fileName = [filesep, sprintf('Energy%g',energy), filesep, 'carbon_Fluence.txt'];
    fileID = fopen(fullfile(workingDir,folderName, fileName),'w+');
    basicFile = fileread([folderName filesep 'ion_Fluence_basic.txt']);
    fprintf(fileID,'%s\n\n',basicFile);
    fprintf(fileID, ' \n');

    fprintf(fileID, 'd:Sim/HL = %f mm \n', currParams.HL );
    fprintf(fileID, 'd:Sim/Range = %f mm \n', currParams.range);
    fprintf(fileID, 'i:Sim/ZCoarse/ZBins = %i \n', currParams.HL );
    fprintf(fileID, ' \n');
    
    fprintf(fileID,'s:So/MySource/BeamEnergySpectrumType = "Continuous"\n')
    fprintf(fileID, 'dv:So/MySource/BeamEnergySpectrumValues = %d ', length(currParams.Ebins));
    fprintf(fileID, '%.6f ', currParams.Ebins);
    fprintf(fileID, 'MeV\n');
    fprintf(fileID, 'uv:So/MySource/BeamEnergySpectrumWeights= %d ', length(currParams.Espec));
    fprintf(fileID, '%.6f ', currParams.Espec);
    fprintf(fileID, '\n');
    fprintf(fileID, 'd:Sim/MySource/SigmaX = %f mm \n', currParams.sigma);
    fprintf(fileID, 'u:Sim/MySource/SigmaXprime = %f \n', currParams.sigmaPrime);
    fprintf(fileID, 'u:Sim/MySource/CorrelationX = %f \n', currParams.correlation);
    fprintf(fileID, 'i:Sim/MySource/NumberOfHistories = %i \n', nPrimary)
    fprintf(fileID, '\n');
    fprintf(fileID, 'd:Sim/EMax_ion = %f MeV \n',currParams.Emax );
    fprintf(fileID, 'd:Sim/EMin_ion = %f MeV \n', 0);
    fprintf(fileID, 'i:Sim/Ebins_ion = %i \n',currParams.Ebins_scoring);

    fclose(fileID);

end

function writeCallFile(workingDir, folderName, energy)
    folderName = [folderName, filesep, sprintf('Energy%g',energy)]; 
    fileName = '\submit_runTOPASCluster';
    fileID = fopen(fullfile(workingDir,folderName, fileName),'w+');
    fprintf(fileID, '#!/bin/bash \nsbatch <<EOT  \n#!/bin/bash \n');
    fprintf(fileID, '#SBATCH --job-name=topas_Seckler_energy%g_run$1',energy );
    fprintf(fileID, '#SBATCH --nodes=1 \n#SBATCH --ntasks=1  \n#SBATCH --cpus-per-task=8 \n');
    fprintf(fileID, '#SBATCH --output=topas_Seckler_energy%g_run$1\n',energy );
    fprintf(fileID, 'pwd; hostname; date \n \n');
    fprintf(fileID, 'export QT_QPA_PLATFORM_PLUGIN_PATH=/opt/topas/OpenTOPAS-main/Frameworks \n');
    fprintf(fileID, 'export TOPAS_G4_DATA_DIR=/opt/GEANT4_11.1.3/G4DATA \n');
    fprintf(fileID, 'export LD_LIBRARY_PATH=/opt/topas/OpenTOPAS-main/lib64:$LD_LIBRARY_PATH \n');
    fprintf(fileID, 'export LD_LIBRARY_PATH=/opt/GEANT4_11.1.3/geant4-install/lib64:$LD_LIBRARY_PATH \n');
    fprintf(fileID, 'export TOPAS_G4_DATA_DIR=/opt/GEANT4_11.1.3/G4DATA \n');
    fprintf(fileID, '/opt/topas/OpenTOPAS-main/bin/topas carbon_DoseAndRBE.txt > carbon_DoseAndRBE.out 2> carbon_DoseAndRBE.err \n');
    fprintf(fileID, '/opt/topas/OpenTOPAS-main/bin/topas carbon_Fluence.txt > carbon_Fluence.out 2> carbon_Fluence.err \n');
    fprintf(fileID, '\nexit 0 \nEOT\n');
    fclose(fileID);
end



