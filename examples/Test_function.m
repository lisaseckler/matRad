%% Testing the function
clear all; close all; clc;

%% Proton Optimization
matRad_rc;

load("PROSTATE.mat")

% cst{1,6}{1} = struct(DirtyDoseObjectives.matRad_SquaredDeviationDirtyDose(400,0));
% 
% cst{16,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredDeviationDirtyDose(400,0));
% 
% cst{15,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose (100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 90;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
% dij = matRad_calcDirtyDose(4,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

index = [79 69 127];
LET_thres = 4;
add = [1 0 0];
% dirtyTest = matRad_plotLETbeamletSpectrumInVoxel(index, ct, 1, dij, resultGUI,50,1);
% dirtyTest = matRad_returnDirtyandCleanDose(index,ct,1,dij,resultGUI,LET_thres,1,[],[]);
dirtyTest3 = matRad_DoseComparison(20,index,ct,1,dij,resultGUI,add,LET_thres,1,1);


%% User Interface
matRadGUI

%% Two Beams
%% Proton Optimization
matRad_rc;

load('BOXPHANTOM_TINY (1).mat')

% matRad_calcDirtyDose with higher resolution is not working!
% ct.resolution.x = 6;
% ct.resolution.y = 6;
% ct.resolution.z = 6;

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 180];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'dirtyDose'; 
modelName     = 'constant'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation

dij = matRad_calcParticleDose(ct,stf,pln,cst);
% dij.ctGrid.resolution.x = 6;
% dij.ctGrid.resolution.y = 6;
% dij.ctGrid.resolution.z = 6;

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(4,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

%% Second dose calculation
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

%% Function Test
index = [75 75 72];
LET_thres = 2.5;
add = [0 1 0];
dirtyTest = matRad_plotLETbeamletSpectrumInVoxel(index, ct, 1, dij, resultGUI,50,1);
dirtyTest2 = matRad_returnDirtyandCleanDose(index,ct,1,dij,resultGUI,LET_thres,1,[],[]);
dirtyTest3 = matRad_DoseComparison(20,index,ct,1,dij,resultGUI,add,LET_thres,2,1);

%% Function loop
index = [74 75 72];
thres = [80 120 80];
add = [0 1 0];
for i =1: 80
    index(2)= index(2) +1;
      letspectrum(i) =matRad_plotLETbeamletSpectrumInVoxel(thres,index,ct,1,dij,resultGUI,add,[],4,6,[],1);

end

%% loop

index = input("Please define the index of your new voxel. Use this: [x y z]. index = ");
bins = input("If you want to use a certain number of bins for your histogram, please define it here: bins = ");
LET_thres = input("If you want to set a threshold to compute dirty and clean dose, please indicate one here: LET_thres = ");
maxdirtydose = input("If you want to set a maximum LET up to which the histogram should be plotted, please indicate one here: maxdirtydose = ");
displayfigures = input("If you want to see the plots, specify displayfigures as 1. displayfigures = ");

while index ~ [];
    matRad_plotLETbeamletSpectrumInVoxel(index,ct,1,dij,resultGUI,bins,LET_thres,maxdirtydose,displayfigures);
    index = input("Please define the index of your new voxel. Use this: [x y z]. index = ");
end

%% Test
index = [80 40 80];
add = [0 1 0];
doseComparison = matRad_DoseComparison(80,index,ct,1,dij,resultGUI,add,4,2,1);


%% Carbon Ions
matRad_rc;
load('PROSTATE.mat')

pln.radiationMode = 'carbon';            
pln.machine       = 'Generic';

modelName           = 'LEM';
quantityOpt         = 'RBExD';   

% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 180];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario                                            

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%Let's also calculate the LET
pln.propDoseCalc.calcLET = true;

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Helium

matRad_rc
load('BOXPHANTOM_TINY (1).mat');

pln.radiationMode = 'helium';        
pln.machine       = 'Generic';

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 90;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

quantityOpt   = 'RBExD';            % either  physicalDose / effect / RBExD
modelName     = 'HEL'; 

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');

stf = matRad_generateStf(ct,cst,pln);

dij = matRad_calcParticleDose(ct,stf,pln,cst);

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

