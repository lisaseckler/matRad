
load BeamOptics.mat
load EnergySpectrumData.mat


%% Simulation Loop

fileName = 'run_all';
fileID = fopen( fileName,'w+');
for i = 1:numel(E)
    
   
    fprintf(fileID, 'cd Energy%g \n', E(i));
    fprintf(fileID, 'bash submit_runTOPASCluster \n')
    fprintf(fileID, 'cd .. \n')
  
end
   
disp('All Simulations Files Written')
fclose all;

%%