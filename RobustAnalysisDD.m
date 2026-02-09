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
cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));
cst{4,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

%% set pln
% only for protons, then only for photons and combine the plan in the end
% meta information for treatment plan (1) 
pln_DD1(1).numOfFractions  = 5;
pln_DD1(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_DD1(1).machine         = 'Generic';

% beam geometry settings
pln_DD1(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_DD1(1).propStf.gantryAngles    = 0;%[-45 45];%[45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_DD1(1).propStf.couchAngles     = zeros(numel(pln_DD1(1).propStf.gantryAngles),1); % [?] ; 
pln_DD1(1).propStf.numOfBeams      = numel(pln_DD1(1).propStf.gantryAngles);
pln_DD1(1).propStf.isoCenter       = ones(pln_DD1(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln_DD1(1).propDoseCalc.calcLET = 1;

pln_DD1(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_DD1(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_DD1(1).propOpt.spatioTemp      = 0;
pln_DD1(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_DD1(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_DD1(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_DD1(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

% worst case scenario for robustness analysis
scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve scenarios for dose calculation and optimziation
pln_DD1(1).multScen = matRad_multScen(ct,scenGenType);

% retrieve bio model parameters
pln_DD1(1).bioParam = matRad_bioModel(pln_DD1(1).radiationMode,quantityOpt, modelName);
%%
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);
%% stf generation single

stf_LxD100_1beam_protons = matRad_generateStf(ct,cst,pln_1beam_DD100(1));

%% multiple scenarios
modality = 'protons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_LxD100_2beams';
multScen_LxD2beams_oldStf = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_2beams_L100(1:2),pln_LxD2(1),cst,modality,[], saveDirectory);

%%
% meta information for treatment plan (2) 
pln_DD1(2).numOfFractions  = 25;
pln_DD1(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_DD1(2).machine         = 'Generic';

% beam geometry settings
pln_DD1(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_DD1(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_DD1(2).propStf.couchAngles     = zeros(numel(pln_DD1(2).propStf.gantryAngles),1);  % [?] ; 
pln_DD1(2).propStf.numOfBeams      = numel(pln_DD1(2).propStf.gantryAngles);
pln_DD1(2).propStf.isoCenter       = ones(pln_DD1(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln_DD1(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_DD1(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_DD1(2).propOpt.spatioTemp      = 0;
pln_DD1(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_DD1(2).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_DD1(2).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_DD1(2).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_DD1(2).bioParam = matRad_bioModel(pln_DD1(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_DD1(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
%%
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);
%% stf generation single

stf_DD100_0_photons = matRad_generateStf(ct,cst,pln_DD1(2));

%% multiple scenarios
modality = 'photons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Photons_DD_100_0';
multScen_LxDphotons = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_1beam_L100(2:end),pln_DD1(2),cst,modality,[], saveDirectory);


%% get meta information
% Photons 9 beams
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Photons_DD_100_0';
scenIdxToInclude = 1:9;
meta_photons = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons DD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_DD100_1beam';
scenIdxToInclude = 1:9;
meta_DD1 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons DD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_DD100_2beam';
scenIdxToInclude = 1:9;
meta_DD2 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons DD -45° and 0° and 45°
% saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_3beams';
% scenIdxToInclude = 1:9;
% meta_DD3 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons LxD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_LxD100_1beam';
scenIdxToInclude = 1:9;
meta_LxD1 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons LxD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_LxD100_2beams';
scenIdxToInclude = 1:9;
meta_LxD2_new = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons LxD -45° and 0° and 45°
% saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_3beams';
% scenIdxToInclude = 1:9;
% meta_LxD3 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

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
    physD_photons{scen} = matRad_loadScenariosFromMeta(saveDir, meta_photons(scen), dij_photons, 0); 
end

% protons DD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_DD100_1beam';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_DD100_1beam\dijTemplate.mat")
dij_protonsDD1 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physD_protonsDD1{scen} = matRad_loadScenariosFromMeta(saveDir, meta_DD1(scen), dij_protonsDD1, 0); 
end

% protons DD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_DD100_2beam';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_DD100_2beam\dijTemplate.mat")
dij_protonsDD2 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physD_protonsDD2{scen} = matRad_loadScenariosFromMeta(saveDir, meta_DD2(scen), dij_protonsDD2, 0); 
end

% protons DD -45° and 0° and 45°
% saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_3beams';
% load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_DD100_3beams\dijTemplate.mat")
% dij_protonsDD3 = dijTemplate;
% stringLength = 0;
% nScenarios = size(scenIdxToInclude,2);
% for scen = scenIdxToInclude
%     fprintf(repmat('\b',1,stringLength));
%     stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
%     physicalDoseDij_protonsDD3{scen} = matRad_loadScenariosFromMeta(saveDir, meta_DD3(scen), dij_protonsDD3, 0); 
% end

% protons LxD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_LxD100_1beam';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_LxD100_1beam\dijTemplate.mat")
dij_protonsLxD1 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physD_protonsLxD1{scen} = matRad_loadScenariosFromMeta(saveDir, meta_LxD1(scen), dij_protonsLxD1, 0); 
end

% protons LxD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_LxD100_2beams';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_oldStf_LxD100_2beams\dijTemplate.mat")
dij_protonsLxD2_new = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    physD_protonsLxD2_new{scen} = matRad_loadScenariosFromMeta(saveDir, meta_LxD2_new(scen), dij_protonsLxD2_new, 0); 
end

% protons LxD -45° and 0° and 45°
% saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_3beams';
% load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD\Protons_LxD100_3beams\dijTemplate.mat")
% dij_protonsLxD3 = dijTemplate;
% stringLength = 0;
% nScenarios = size(scenIdxToInclude,2);
% for scen = scenIdxToInclude
%     fprintf(repmat('\b',1,stringLength));
%     stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
%     physicalDoseDij_protonsLxD3{scen} = matRad_loadScenariosFromMeta(saveDir, meta_LxD3(scen), dij_protonsLxD3, 0); 
% end

%% Weighting physicalDose

for scen = scenIdxToInclude
    photonsPhysD{scen} = physD_photons{1,scen}.physicalDose{1,1} * result_overDD_1beam{1,2}.w;
    protonDD1PhysD{scen} = physD_protonsDD1{1,scen}.physicalDose{1,1} * result_overDD_1beam{1,1}.w;
    protonDD2PhysD{scen} = physD_protonsDD2{1,scen}.physicalDose{1,1} * result_overDD_2beams{1,1}.w;
    protonLxD1PhysD{scen} = physD_protonsLxD1{1,scen}.physicalDose{1,1} * result_overLxD100_1beam{1,1}.w;
    protonLxD2PhysD{scen} = physD_protonsLxD2_new{1,scen}.physicalDose{1,1} * result_overLxD100_2beams{1,1}.w;
end

%%
% Plan Wrapper
plnJO_1beam_DD100 = matRad_plnWrapper(pln_DD1);
% Stf Wrapper
stf_1beam_DD100 = matRad_stfWrapper(ct,cst,pln_DD1);

% Dij Calculation
dij_1beam_DD100 = matRad_calcCombiDose(ct,stf_1beamDD100_test,plnJO_1beam_DD100,cst,false);
% Dirty Dose Calculation
dij_1beam_DD100 = matRad_calcDirtyDose(2,dij_1beam_DD100,pln_DD1);
dij_1beam_DD100 = matRad_calcmLETDose(dij_1beam_DD100,pln_DD1);
dij_1beam_DD100.precon = 0;
%%
[result_overLxD60_1beam,optimizer_overLxD60_1beam] = matRad_fluenceOptimizationJO(dij_1beam_DD100,cst,plnJO_1beam_DD100);

%% Reshaping

for scen = scenIdxToInclude
    PhotonPhysD{scen} = reshape(photonsPhysD{1}, 167, 167, 107);
    ProtonDD1PhysD{scen} = reshape(protonDD1PhysD{scen}, 167, 167, 107);
    ProtonDD2PhysD{scen} = reshape(protonDD2PhysD{scen}, 167, 167, 107);
    ProtonLxD1PhysD{scen} = reshape(protonLxD1PhysD{scen}, 167, 167, 107);
    ProtonLxD2PhysD{scen} = reshape(protonLxD2PhysD{scen}, 167, 167, 107);

end

%% Combining
for scen = scenIdxToInclude
    CombiDD1_sameScen{scen} = 5 * ProtonDD1PhysD{scen} + 25 * PhotonPhysD{scen};
    CombiDD2_sameScen{scen} = 5 * ProtonDD2PhysD{scen} + 25 * PhotonPhysD{scen};
    CombiLxD1_sameScen{scen} = 5 * ProtonLxD1PhysD{scen} + 25 * PhotonPhysD{scen};
    CombiLxD2_sameScen{scen} = 5 * ProtonLxD2PhysD{scen} + 25 * PhotonPhysD{scen};
end
 

%% DVH
for scen = scenIdxToInclude
    DVH_DD1{scen} = matRad_calcDVH(cst,CombiDD1_sameScen{scen});
    DVH_DD2{scen} = matRad_calcDVH(cst,CombiDD2_sameScen{scen});
    DVH_LxD1{scen} = matRad_calcDVH(cst,CombiLxD1_sameScen{scen});
    DVH_LxD2{scen} = matRad_calcDVH(cst,CombiLxD2_sameScen{scen});
end

%% Plotting DVH
color = ['b', 'm', 'g' , 'r'];

figure
structure = [1 2 3];
for i = [1 2 3]
    plot(DVH_DD1{1}(structure(i)).doseGrid,DVH_DD1{1}(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',2)
    hold on
    plot(DVH_DD2{1}(structure(i)).doseGrid,DVH_DD2{1}(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',2)
    hold on
    plot(DVH_LxD1{1}(structure(i)).doseGrid,DVH_LxD1{1}(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',2)
    hold on
    plot(DVH_LxD2{1}(structure(i)).doseGrid,DVH_LxD2{1}(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',2)
    hold on    
end
