%% Robustness analysis for Paper 
% calculated with DirtyDoseComparison.m

matRad_rc;
clear
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 500000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-05;
load TG119.mat

%% add margin
cube = zeros(ct.cubeDim);
cube(cst{1,4}{1}) = 1;
vResolution = ct.resolution;
vMargin = [];
vMargin.x = 5;
vMargin.y = 5;
vMargin.z = 5;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin,1);

cst{4,1}    = 3;
cst{4,2}    = 'Core_Big';
cst{4,3}    = 'OAR';
cst{4,4}{1} = find(mVOIEnlarged);
cst{4,5}    = cst{1,5};

%% primary objectives
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,60));
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));
cst{3,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));

%% LMO objectives
%cst{2,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,30));
cst{1,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,0));
cst{4,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,0));

%% set pln
% only for protons, then only for photons and combine the plan in the end
% meta information for treatment plan (1) 
pln_3beams_DD100(1).numOfFractions  = 5;
pln_3beams_DD100(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_3beams_DD100(1).machine         = 'Generic';

% beam geometry settings
pln_3beams_DD100(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_3beams_DD100(1).propStf.gantryAngles    = 0; %[-45 45];%[45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_3beams_DD100(1).propStf.couchAngles     = zeros(numel(pln_3beams_DD100(1).propStf.gantryAngles),1); % [?] ; 
pln_3beams_DD100(1).propStf.numOfBeams      = numel(pln_3beams_DD100(1).propStf.gantryAngles);
pln_3beams_DD100(1).propStf.isoCenter       = ones(pln_3beams_DD100(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln_3beams_DD100(1).propDoseCalc.calcLET = 1;

pln_3beams_DD100(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_3beams_DD100(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_3beams_DD100(1).propOpt.spatioTemp      = 0;
pln_3beams_DD100(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_3beams_DD100(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_3beams_DD100(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_3beams_DD100(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

% worst case scenario for robustness analysis
scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve scenarios for dose calculation and optimziation
pln_3beams_DD100(1).multScen = matRad_multScen(ct,scenGenType);

% retrieve bio model parameters
pln_3beams_DD100(1).bioParam = matRad_bioModel(pln_3beams_DD100(1).radiationMode,quantityOpt, modelName);
%%
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);
%% stf generation single

stf_LxD100_1beam_protons = matRad_generateStf(ct,cst,pln_1beam_DD100(1));

%% multiple scenarios
modality = 'protons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_1beam';
multScen_LxD1beam = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_LxD100_1beam_protons,pln_1beam_DD100(1),cst,modality,[], saveDirectory);

%%
% meta information for treatment plan (2) 
pln_3beams_DD100(2).numOfFractions  = 25;
pln_3beams_DD100(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_3beams_DD100(2).machine         = 'Generic';

% beam geometry settings
pln_3beams_DD100(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_3beams_DD100(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_3beams_DD100(2).propStf.couchAngles     = zeros(numel(pln_3beams_DD100(2).propStf.gantryAngles),1);  % [?] ; 
pln_3beams_DD100(2).propStf.numOfBeams      = numel(pln_3beams_DD100(2).propStf.gantryAngles);
pln_3beams_DD100(2).propStf.isoCenter       = ones(pln_3beams_DD100(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln_3beams_DD100(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_3beams_DD100(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_3beams_DD100(2).propOpt.spatioTemp      = 0;
pln_3beams_DD100(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_3beams_DD100(2).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_3beams_DD100(2).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_3beams_DD100(2).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_3beams_DD100(2).bioParam = matRad_bioModel(pln_3beams_DD100(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_3beams_DD100(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
%%
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);
%% stf generation single

stf_DD100_0_photons = matRad_generateStf(ct,cst,pln_3beams_DD100(2));

%% multiple scenarios
modality = 'photons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Photons_DD_100_0';
multScen = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_DD100_0_photons,pln_3beams_DD100(2),cst,modality,[], saveDirectory);


%% get meta information
% Photons 9 beams
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Photons_DD_100_0';
scenIdxToInclude = 1:9;
meta_photons = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons DD 0°
% saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_0';
% scenIdxToInclude = 1:9;
% meta_DD1 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons DD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_2beams';
scenIdxToInclude = 1:9;
meta_DD2 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons DD -45° and 0° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_3beams';
scenIdxToInclude = 1:9;
meta_DD3 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons LxD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_1beam';
scenIdxToInclude = 1:9;
meta_LxD1 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons LxD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_2beams';
scenIdxToInclude = 1:9;
meta_LxD2 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons LxD -45° and 0° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_3beams';
scenIdxToInclude = 1:9;
meta_LxD3 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

%% load meta information

% photons 9 beams
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Photons_DD_100_0';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Photons_DD_100_0\dijTemplate.mat")
dij_photons = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physicalDoseDij_photons{scen} = matRad_loadScenariosFromMeta(saveDir, meta_photons(scen), dij_photons, 0); 
end

% protons DD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_0';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_0\dijTemplate.mat")
dij_protonsDD1 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physicalDoseDij_protonsDD1{scen} = matRad_loadScenariosFromMeta(saveDir, meta_DD1(scen), dij_protonsDD1, 0); 
end

% protons DD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_2beams';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_2beams\dijTemplate.mat")
dij_protonsDD2 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physicalDoseDij_protonsDD2{scen} = matRad_loadScenariosFromMeta(saveDir, meta_DD2(scen), dij_protonsDD2, 0); 
end

% protons DD -45° and 0° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_3beams';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_3beams\dijTemplate.mat")
dij_protonsDD3 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physicalDoseDij_protonsDD3{scen} = matRad_loadScenariosFromMeta(saveDir, meta_DD3(scen), dij_protonsDD3, 0); 
end

% protons LxD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_1beam';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_1beam\dijTemplate.mat")
dij_protonsLxD1 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physicalDoseDij_protonsLxD1{scen} = matRad_loadScenariosFromMeta(saveDir, meta_LxD1(scen), dij_protonsLxD1, 0); 
end

% protons LxD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_2beams';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_2beams\dijTemplate.mat")
dij_protonsLxD2 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physicalDoseDij_protonsLxD2{scen} = matRad_loadScenariosFromMeta(saveDir, meta_LxD2(scen), dij_protonsLxD2, 0); 
end

% protons LxD -45° and 0° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_3beams';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_3beams\dijTemplate.mat")
dij_protonsLxD3 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physicalDoseDij_protonsLxD3{scen} = matRad_loadScenariosFromMeta(saveDir, meta_LxD3(scen), dij_protonsLxD3, 0); 
end


%%
% Plan Wrapper
plnJO_1beam_DD100 = matRad_plnWrapper(pln_3beams_DD100);
% Stf Wrapper
stf_1beam_DD100 = matRad_stfWrapper(ct,cst,pln_3beams_DD100);

% Dij Calculation
dij_1beam_DD100 = matRad_calcCombiDose(ct,stf_1beamDD100_test,plnJO_1beam_DD100,cst,false);
% Dirty Dose Calculation
dij_1beam_DD100 = matRad_calcDirtyDose(2,dij_1beam_DD100,pln_3beams_DD100);
dij_1beam_DD100 = matRad_calcmLETDose(dij_1beam_DD100,pln_3beams_DD100);
dij_1beam_DD100.precon = 0;
%%
[result_overLxD60_1beam,optimizer_overLxD60_1beam] = matRad_fluenceOptimizationJO(dij_1beam_DD100,cst,plnJO_1beam_DD100);
