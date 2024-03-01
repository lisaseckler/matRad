%% Testing the function
clear; close all; clc;

%% Proton Optimization 2 beam without dirty Dose
matRad_rc;
clear("dij","pln","resultGUI","ct","cst","stf")
load("PROSTATE.mat")

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIRef = resultGUI;

save("Reference_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 100
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver100 = resultGUI;

save("BladderOver100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 300
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver300 = resultGUI;

save("BladderOver300_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 500
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver500 = resultGUI;

save("BladderOver500_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 700
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(700,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver700 = resultGUI;

save("BladderOver700_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 900
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver900 = resultGUI;

save("BladderOver900_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 100 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder100 = resultGUI;

save("TargetUnder100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 300 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder300 = resultGUI;

save("TargetUnder300_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 500 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder500 = resultGUI;

save("TargetUnder500_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 700 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(700,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder700 = resultGUI;

save("TargetUnder700_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 900 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(900,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder900 = resultGUI;

save("TargetUnder900_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 100 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(100,20));
cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(10,10));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 120 240 270];
pln.propStf.couchAngles   = [0 0 0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultTarget = resultGUI;

% save("BodyMean100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')
% 
% totalphysDose = sum(resultGUIMean100.physicalDose(:));
% physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean100.physicalDose);
% Mean100 = matRad_calcQualityIndicators(cst,pln,resultGUIMean100.dirtyDose);
% sumDirtyDose = sum(resultGUIMean100.dirtyDose(cst{8,4}{1}));
% maxDirtyDose = max(resultGUIMean100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
% save("BodyMean100_2BeamConst.mat","totalphysDose","physDose","Mean100","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 100 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(100,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean100 = resultGUI;

save("BodyMean100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean100.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean100.physicalDose);
Mean100 = matRad_calcQualityIndicators(cst,pln,resultGUIMean100.dirtyDose);
sumDirtyDose = sum(resultGUIMean100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean100_2BeamConst.mat","totalphysDose","physDose","Mean100","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 300 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(300,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean300 = resultGUI;

save("BodyMean300_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean300.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean300.physicalDose);
Mean300 = matRad_calcQualityIndicators(cst,pln,resultGUIMean300.dirtyDose);
sumDirtyDose = sum(resultGUIMean300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean300_2BeamConst.mat","totalphysDose","physDose","Mean300","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 500 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(500,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean500 = resultGUI;

save("BodyMean500_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean500.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean500.physicalDose);
Mean500 = matRad_calcQualityIndicators(cst,pln,resultGUIMean500.dirtyDose);
sumDirtyDose = sum(resultGUIMean500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean500_2BeamConst.mat","totalphysDose","physDose","Mean500","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 700 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(700,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean700 = resultGUI;

save("BodyMean700_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean700.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean700.physicalDose);
Mean700 = matRad_calcQualityIndicators(cst,pln,resultGUIMean700.dirtyDose);
sumDirtyDose = sum(resultGUIMean700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean700_2BeamConst.mat","totalphysDose","physDose","Mean700","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 900 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(900,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean900 = resultGUI;

save("BodyMean900_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean900.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean900.physicalDose);
Mean900 = matRad_calcQualityIndicators(cst,pln,resultGUIMean900.dirtyDose);
sumDirtyDose = sum(resultGUIMean900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean900_2BeamConst.mat","totalphysDose","physDose","Mean900","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder Overdosing penalty 100 Target underdosing dirty Dose with penalty 100 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver100Under100 = resultGUI;

save("Over100Under100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean900.physicalDose);
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean900.physicalDose);
Mean900 = matRad_calcQualityIndicators(cst,pln,resultGUIMean900.dirtyDose);
sumDirtyDose = sum(resultGUIMean900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean900_2BeamConst.mat","totalphysDose","physDose","Mean900","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder Overdosing penalty 100 Target underdosing dirty Dose with penalty 300 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver100Under300 = resultGUI;

save("Over100Under300_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 100 Target underdosing dirty Dose with penalty 500 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver100Under500 = resultGUI;

save("Over100Under500_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 100 Target underdosing dirty Dose with penalty 700 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(700,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver100Under700 = resultGUI;

save("Over100Under700_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 100 Target underdosing dirty Dose with penalty 900 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(900,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver100Under900 = resultGUI;

save("Over100Under900_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 300 Target underdosing dirty Dose with penalty 100 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver300Under100 = resultGUI;

save("Over300Under100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 300 Target underdosing dirty Dose with penalty 300 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver300Under300 = resultGUI;

save("Over300Under300_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 300 Target underdosing dirty Dose with penalty 500 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver300Under500 = resultGUI;

save("Over300Under500_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 300 Target underdosing dirty Dose with penalty 700 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(700,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver300Under700 = resultGUI;

save("Over300Under700_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 300 Target underdosing dirty Dose with penalty 900 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(900,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver300Under900 = resultGUI;

save("Over300Under900_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 500 Target underdosing dirty Dose with penalty 100 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver500Under100 = resultGUI;

save("Over500Under100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 500 Target underdosing dirty Dose with penalty 300 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver500Under300 = resultGUI;

save("Over500Under300_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 500 Target underdosing dirty Dose with penalty 500 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver500Under500 = resultGUI;

save("Over500Under500_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 500 Target underdosing dirty Dose with penalty 700 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(700,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver500Under700 = resultGUI;

save("Over500Under700_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 500 Target underdosing dirty Dose with penalty 900 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver500Under900 = resultGUI;

save("Over500Under900_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 700 Target underdosing dirty Dose with penalty 100 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(700,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver700Under100 = resultGUI;

save("Over700Under100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 700 Target underdosing dirty Dose with penalty 300 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(700,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver700Under300 = resultGUI;

save("Over700Under300_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 700 Target underdosing dirty Dose with penalty 500 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(700,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver700Under500 = resultGUI;

save("Over700Under500_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 700 Target underdosing dirty Dose with penalty 700 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(700,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(700,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver700Under700 = resultGUI;

save("Over700Under700_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 700 Target underdosing dirty Dose with penalty 900 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(900,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(700,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver700Under900 = resultGUI;

save("Over700Under900_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 900 Target underdosing dirty Dose with penalty 100 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver900Under100 = resultGUI;

save("Over900Under100_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 900 Target underdosing dirty Dose with penalty 300 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver900Under300 = resultGUI;

save("Over900Under300_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 900 Target underdosing dirty Dose with penalty 500 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver900Under500 = resultGUI;

save("Over900Under500_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 900 Target underdosing dirty Dose with penalty 700 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(700,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver900Under700 = resultGUI;

save("Over900Under700_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')

%% Proton Optimization 2 beam Bladder Overdosing penalty 900 Target underdosing dirty Dose with penalty 900 dmax 60
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(900,60));
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver900Under900 = resultGUI;

save("Over900Under900_2BeamConst.mat","resultGUI","pln","dij","ct","cst",'-mat')
%% Ref Plan DVH
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIRef.physicalDose);
Ref = matRad_calcQualityIndicators(cst,pln,resultGUIRef.dirtyDose);

doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
sumDirtyDose = sum(resultGUIRef.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIRef.dirtyDose(cst{8,4}{1})); 
totalphysDose = sum(resultGUIRef.physicalDose(:));
save("Reference_2BeamConst.mat","totalphysDose","physDose","Ref","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);

plane = 3;
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIRef.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
% (axesHandle,ct,cst,cubeIdx,dose,plane,slice,thresh,alpha,contourColorMap,doseColorMap,doseWindow,doseIsoLevels,voiSelection,colorBarLabel,boolPlotLegend,varargin)
title('original plan')

%% Beam 2 Over 100 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver100.physicalDose);
Over100 = matRad_calcQualityIndicators(cst,pln,resultGUIOver100.dirtyDose);
sumDirtyDose = sum(resultGUIOver100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIOver100.physicalDose(:));
save("BladderOver100_2BeamConst.mat","totalphysDose","physDose","Over100","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver100.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver100.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 300 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver300.physicalDose);
Over300 = matRad_calcQualityIndicators(cst,pln,resultGUIOver300.dirtyDose);
sumDirtyDose = sum(resultGUIOver300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIOver300.physicalDose(:));
save("BladderOver300_2BeamConst.mat","totalphysDose","physDose","Over300","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver300.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver300.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 500 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver500.physicalDose);
Over500 = matRad_calcQualityIndicators(cst,pln,resultGUIOver500.dirtyDose);
sumDirtyDose = sum(resultGUIOver500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIOver500.physicalDose(:));
save("BladderOver500_2BeamConst.mat","totalphysDose","physDose","Over500","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver500.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver500.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 700 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver700.physicalDose);
Over700 = matRad_calcQualityIndicators(cst,pln,resultGUIOver700.dirtyDose);
sumDirtyDose = sum(resultGUIOver700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIOver700.physicalDose(:));
save("BladderOver700_2BeamConst.mat","totalphysDose","physDose","Over700","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver700.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver700.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 900 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver900.physicalDose);
Over900 = matRad_calcQualityIndicators(cst,pln,resultGUIOver900.dirtyDose);
sumDirtyDose = sum(resultGUIOver900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIOver900.physicalDose(:));
save("BladderOver900_2BeamConst.mat","totalphysDose","physDose","Over900","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver900.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver900.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Under 100 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder100.physicalDose);
Under100 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder100.dirtyDose);
sumDirtyDose = sum(resultGUIUnder100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIUnder100.physicalDose(:));
save("TargetUnder100_2BeamConst.mat","totalphysDose","physDose","Under100","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = 45;

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIUnder100.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIUnder100.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Under 300 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder300.physicalDose);
Under300 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder300.dirtyDose);
sumDirtyDose = sum(resultGUIUnder300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIUnder300.physicalDose(:));
save("TargetUnder300_2BeamConst.mat","totalphysDose","physDose","Under300","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = 45;

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIUnder300.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIUnder300.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Under 500 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder500.physicalDose);
Under500 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder500.dirtyDose);
sumDirtyDose = sum(resultGUIUnder500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIUnder500.physicalDose(:));
save("TargetUnder500_2BeamConst.mat","totalphysDose","physDose","Under500","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = 45;

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIUnder500.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIUnder500.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Under 700 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder700.physicalDose);
Under700 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder700.dirtyDose);
sumDirtyDose = sum(resultGUIUnder700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIUnder700.physicalDose(:));
save("TargetUnder700_2BeamConst.mat","totalphysDose","physDose","Under700","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = 45;

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIUnder700.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIUnder700.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Under 900 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder900.physicalDose);
Under900 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder900.dirtyDose);
sumDirtyDose = sum(resultGUIUnder900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
totalphysDose = sum(resultGUIUnder900.physicalDose(:));
save("TargetUnder900_2BeamConst.mat","totalphysDose","physDose","Under900","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = 45;

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIUnder900.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIUnder900.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 100 Under 100 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under100.physicalDose);
Over100Under100 = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under100.dirtyDose);
sumDirtyDose = sum(resultGUIOver100Under100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver100Under100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over100Under100_2BeamConst.mat","physDose","Over100Under100","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver100Under100.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver100Under100.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 100 Under 300 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under300.physicalDose);
Over100Under300 = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under300.dirtyDose);
sumDirtyDose = sum(resultGUIOver100Under300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver100Under300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over100Under300_2BeamConst.mat","physDose","Over100Under300","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver100Under300.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver100Under300.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 100 Under 500 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under500.physicalDose);
Over100Under500 = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under500.dirtyDose);
sumDirtyDose = sum(resultGUIOver100Under500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver100Under500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over100Under500_2BeamConst.mat","physDose","Over100Under500","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver100Under500.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver100Under500.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 100 Under 700 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under700.physicalDose);
Over100Under700 = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under700.dirtyDose);
sumDirtyDose = sum(resultGUIOver100Under700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver100Under700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over100Under700_2BeamConst.mat","physDose","Over100Under700","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver100Under700.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver100Under700.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 100 Under 900 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under900.physicalDose);
Over100Under900 = matRad_calcQualityIndicators(cst,pln,resultGUIOver100Under900.dirtyDose);
sumDirtyDose = sum(resultGUIOver100Under900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver100Under900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over100Under900_2BeamConst.mat","physDose","Over100Under900","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver100Under900.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver100Under900.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 300 Under 100 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under100.physicalDose);
Over300Under100 = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under100.dirtyDose);
sumDirtyDose = sum(resultGUIOver300Under100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver300Under100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over300Under100_2BeamConst.mat","physDose","Over300Under100","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver300Under100.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver300Under100.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 300 Under 300 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under300.physicalDose);
Over300Under300 = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under300.dirtyDose);
sumDirtyDose = sum(resultGUIOver300Under300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver300Under300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over300Under300_2BeamConst.mat","physDose","Over300Under300","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver300Under300.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver300Under300.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 300 Under 500 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under500.physicalDose);
Over300Under500 = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under500.dirtyDose);
sumDirtyDose = sum(resultGUIOver300Under500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver300Under500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over300Under500_2BeamConst.mat","physDose","Over300Under500","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver300Under500.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver300Under500.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 300 Under 700 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under700.physicalDose);
Over300Under700 = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under700.dirtyDose);
sumDirtyDose = sum(resultGUIOver300Under700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver300Under700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over300Under700_2BeamConst.mat","physDose","Over300Under700","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver300Under700.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver300Under700.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 300 Under 900 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under900.physicalDose);
Over300Under900 = matRad_calcQualityIndicators(cst,pln,resultGUIOver300Under900.dirtyDose);
sumDirtyDose = sum(resultGUIOver300Under900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver300Under900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over300Under900_2BeamConst.mat","physDose","Over300Under900","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver300Under900.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver300Under900.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 500 Under 100 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under100.physicalDose);
Over500Under100 = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under100.dirtyDose);
sumDirtyDose = sum(resultGUIOver500Under100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver500Under100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over500Under100_2BeamConst.mat","physDose","Over500Under100","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver500Under100.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver500Under100.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 500 Under 300 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under300.physicalDose);
Over500Under300 = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under300.dirtyDose);
sumDirtyDose = sum(resultGUIOver500Under300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver500Under300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over500Under300_2BeamConst.mat","physDose","Over500Under300","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver500Under300.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver500Under300.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 500 Under 500 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under500.physicalDose);
Over500Under500 = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under500.dirtyDose);
sumDirtyDose = sum(resultGUIOver500Under500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIRef.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over500Under500_2BeamConst.mat","physDose","Over500Under500","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver500Under500.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver500Under500.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 500 Under 700 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under700.physicalDose);
Over500Under700 = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under700.dirtyDose);
sumDirtyDose = sum(resultGUIOver500Under700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver500Under700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over500Under700_2BeamConst.mat","physDose","Over500Under700","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver500Under700.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver500Under700.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 500 Under 900 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under900.physicalDose);
Over500Under900 = matRad_calcQualityIndicators(cst,pln,resultGUIOver500Under900.dirtyDose);
sumDirtyDose = sum(resultGUIOver500Under900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver500Under900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over500Under900_2BeamConst.mat","physDose","Over500Under900","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver500Under900.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver500Under900.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 700 Under 100 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under100.physicalDose);
Over700Under100 = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under100.dirtyDose);
sumDirtyDose = sum(resultGUIOver700Under100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver700Under100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over700Under100_2BeamConst.mat","physDose","Over700Under100","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver700Under100.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver700Under100.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 700 Under 300 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under300.physicalDose);
Over700Under300 = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under300.dirtyDose);
sumDirtyDose = sum(resultGUIOver700Under300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver700Under300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over700Under300_2BeamConst.mat","physDose","Over700Under300","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver700Under300.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver700Under300.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 700 Under 500 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under500.physicalDose);
Over700Under500 = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under500.dirtyDose);
sumDirtyDose = sum(resultGUIOver700Under500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver700Under500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over700Under500_2BeamConst.mat","physDose","Over700Under500","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver700Under500.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver700Under500.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 700 Under 700 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under700.physicalDose);
Over700Under700 = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under700.dirtyDose);
sumDirtyDose = sum(resultGUIOver700Under700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver700Under700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over700Under700_2BeamConst.mat","physDose","Over700Under700","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver700Under700.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver700Under700.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 700 Under 900 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under900.physicalDose);
Over700Under900 = matRad_calcQualityIndicators(cst,pln,resultGUIOver700Under900.dirtyDose);
sumDirtyDose = sum(resultGUIOver700Under900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver700Under900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over700Under900_2BeamConst.mat","physDose","Over700Under900","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver700Under900.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver700Under900.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 900 Under 100 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under100.physicalDose);
Over900Under100 = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under100.dirtyDose);
sumDirtyDose = sum(resultGUIOver900Under100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver900Under100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over900Under100_2BeamConst.mat","physDose","Over900Under100","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver900Under100.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver900Under100.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 900 Under 300 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under300.physicalDose);
Over900Under300 = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under300.dirtyDose);
sumDirtyDose = sum(resultGUIOver900Under300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver900Under300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over900Under300_2BeamConst.mat","physDose","Over900Under300","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver900Under300.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver900Under300.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 900 Under 500 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under500.physicalDose);
Over900Under500 = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under500.dirtyDose);
sumDirtyDose = sum(resultGUIOver900Under500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver900Under500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over900Under500_2BeamConst.mat","physDose","Over900Under500","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver900Under500.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver900Under500.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 900 Under 700 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under700.physicalDose);
Over900Under700 = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under700.dirtyDose);
sumDirtyDose = sum(resultGUIOver900Under700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver900Under700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over900Under700_2BeamConst.mat","physDose","Over900Under700","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver900Under700.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver900Under700.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')

%% Beam 2 Over 900 Under 900 Comparison
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under900.physicalDose);
Over900Under900 = matRad_calcQualityIndicators(cst,pln,resultGUIOver900Under900.dirtyDose);
sumDirtyDose = sum(resultGUIOver900Under900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver900Under900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Over900Under900_2BeamConst.mat","physDose","Over900Under900","sumDirtyDose","maxDirtyDose",'-mat','-append')

slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);

plane = 3;
doseWindow = [0 max([resultGUIRef.RBExD(:); resultGUIOver100.RBExD(:)])];
diff = resultGUIRef.RBExD-resultGUIOver900Under900.RBExD;

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIOver900Under900.RBExD,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,diff,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')