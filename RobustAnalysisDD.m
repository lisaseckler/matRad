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
cst_DD = cst;
cst_DD{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));
cst_DD{4,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));
%% pln
% % meta information for treatment plan (1) 
pln_DD1_robust(1).numOfFractions  = 5;
pln_DD1_robust(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_DD1_robust(1).machine         = 'Generic';

% beam geometry settings
pln_DD1_robust(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_DD1_robust(1).propStf.gantryAngles    = 0;%[265 285]; %[90 180 270];% [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_DD1_robust(1).propStf.couchAngles     = zeros(numel(pln_DD1_robust(1).propStf.gantryAngles),1); % [?] ; 
pln_DD1_robust(1).propStf.numOfBeams      = numel(pln_DD1_robust(1).propStf.gantryAngles);
pln_DD1_robust(1).propStf.isoCenter       = ones(pln_DD1_robust(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_DD,ct,0);
% optimization settings
pln_DD1_robust(1).propDoseCalc.calcLET = 1;

pln_DD1_robust(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_DD1_robust(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_DD1_robust(1).propOpt.spatioTemp      = 0;
pln_DD1_robust(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_DD1_robust(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_DD1_robust(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_DD1_robust(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
%quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
%=======================================> Model check error in bioModel
%modelName    = 'MCN';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions
modelName = 'constRBE';


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_DD1_robust(1).bioParam = matRad_bioModel(pln_DD1_robust(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_DD1_robust(1).multScen = matRad_multScen(ct,scenGenType);

%pln = pln(1);

% meta information for treatment plan (2) 
pln_DD1_robust(2).numOfFractions  = 25;
pln_DD1_robust(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_DD1_robust(2).machine         = 'Generic';

% beam geometry settings
pln_DD1_robust(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_DD1_robust(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_DD1_robust(2).propStf.couchAngles     = zeros(numel(pln_DD1_robust(2).propStf.gantryAngles),1);  % [?] ; 
pln_DD1_robust(2).propStf.numOfBeams      = numel(pln_DD1_robust(2).propStf.gantryAngles);
pln_DD1_robust(2).propStf.isoCenter       = ones(pln_DD1_robust(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_DD,ct,0);
% optimization settings
pln_DD1_robust(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_DD1_robust(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_DD1_robust(2).propOpt.spatioTemp      = 0;
pln_DD1_robust(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_DD1_robust(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_DD1_robust(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_DD1_robust(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

%quantityOpt  = ['effect'];     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_DD1_robust(2).bioParam = matRad_bioModel(pln_DD1_robust(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_DD1_robust(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst_DD = matRad_prepCst(cst_DD, sparecst);

%pln = pln(2);

% Plan Wrapper
plnJO_DD1 = matRad_plnWrapper(pln_DD1_robust);

%% Stf Wrapper
stf_DD1 = matRad_stfWrapper(ct,cst_DD,plnJO_DD1);

%% Dij Calculation
dij_DD1 = matRad_calcCombiDose(ct,stf_DD1,plnJO_DD1,cst_DD,false);
% Dirty Dose Calculation
dij_DD1 = matRad_calcDirtyDose(2,dij_DD1,pln_DD1_robust);
dij_DD1 = matRad_calcmLETDose(dij_DD1,pln_DD1_robust);

%% compress optimized results into one cube 
%
dij_DD1.precon = 1;
% dij.wInit = [result_pre{1}.w; result_pre{2}.w];
result_DD1 = matRad_fluenceOptimizationJO(dij_DD1,cst_DD,plnJO_DD1);



%% set pln
% only for protons, then only for photons and combine the plan in the end
% meta information for treatment plan (1) 
pln_DD1_robust(1).numOfFractions  = 5;
pln_DD1_robust(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_DD1_robust(1).machine         = 'Generic';

% beam geometry settings
pln_DD1_robust(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_DD1_robust(1).propStf.gantryAngles    = 0;%[-45 45];%[45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_DD1_robust(1).propStf.couchAngles     = zeros(numel(pln_DD1_robust(1).propStf.gantryAngles),1); % [?] ; 
pln_DD1_robust(1).propStf.numOfBeams      = numel(pln_DD1_robust(1).propStf.gantryAngles);
pln_DD1_robust(1).propStf.isoCenter       = ones(pln_DD1_robust(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_DD,ct,0);
% optimization settings
pln_DD1_robust(1).propDoseCalc.calcLET = 1;

pln_DD1_robust(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_DD1_robust(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_DD1_robust(1).propOpt.spatioTemp      = 0;
pln_DD1_robust(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_DD1_robust(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_DD1_robust(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_DD1_robust(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

% worst case scenario for robustness analysis
scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve scenarios for dose calculation and optimziation
pln_DD1_robust(1).multScen = matRad_multScen(ct,scenGenType);

% retrieve bio model parameters
pln_DD1_robust(1).bioParam = matRad_bioModel(pln_DD1_robust(1).radiationMode,quantityOpt, modelName);

plnJO_DD1_robust = matRad_plnWrapper(pln_DD1_robust);

%% multiple scenarios
modality = 'protons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_DD1_new';
multScen_DD1 = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_DD1(1),pln_DD1_robust,cst_DD,modality,[], saveDirectory);

%% pln
% % meta information for treatment plan (1) 
pln_DD2(1).numOfFractions  = 5;
pln_DD2(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_DD2(1).machine         = 'Generic';

% beam geometry settings
pln_DD2(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_DD2(1).propStf.gantryAngles    = [-45 45];%[265 285]; %[90 180 270];% [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_DD2(1).propStf.couchAngles     = zeros(numel(pln_DD2(1).propStf.gantryAngles),1); % [?] ; 
pln_DD2(1).propStf.numOfBeams      = numel(pln_DD2(1).propStf.gantryAngles);
pln_DD2(1).propStf.isoCenter       = ones(pln_DD2(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_DD,ct,0);
% optimization settings
pln_DD2(1).propDoseCalc.calcLET = 1;

pln_DD2(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_DD2(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_DD2(1).propOpt.spatioTemp      = 0;
pln_DD2(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_DD2(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_DD2(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_DD2(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
%quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
%=======================================> Model check error in bioModel
%modelName    = 'MCN';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions
modelName = 'constRBE';


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_DD2(1).bioParam = matRad_bioModel(pln_DD2(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_DD2(1).multScen = matRad_multScen(ct,scenGenType);

%pln = pln(1);

% meta information for treatment plan (2) 
pln_DD2(2).numOfFractions  = 25;
pln_DD2(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_DD2(2).machine         = 'Generic';

% beam geometry settings
pln_DD2(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_DD2(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_DD2(2).propStf.couchAngles     = zeros(numel(pln_DD2(2).propStf.gantryAngles),1);  % [?] ; 
pln_DD2(2).propStf.numOfBeams      = numel(pln_DD2(2).propStf.gantryAngles);
pln_DD2(2).propStf.isoCenter       = ones(pln_DD2(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_DD,ct,0);
% optimization settings
pln_DD2(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_DD2(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_DD2(2).propOpt.spatioTemp      = 0;
pln_DD2(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_DD2(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_DD2(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_DD2(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

%quantityOpt  = ['effect'];     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_DD2(2).bioParam = matRad_bioModel(pln_DD2(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_DD2(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst_DD = matRad_prepCst(cst_DD, sparecst);

%pln = pln(2);

% Plan Wrapper
plnJO_DD2 = matRad_plnWrapper(pln_DD2);

%% Stf Wrapper
stf_DD2 = matRad_stfWrapper(ct,cst_DD,plnJO_DD2);

%% Dij Calculation
dij_DD2 = matRad_calcCombiDose(ct,stf_DD2,plnJO_DD2,cst_DD,false);
% Dirty Dose Calculation
dij_DD2 = matRad_calcDirtyDose(2,dij_DD2,pln_DD2);
dij_DD2 = matRad_calcmLETDose(dij_DD2,pln_DD2);

%% compress optimized results into one cube 
%
dij_DD2.precon = 1;
% dij.wInit = [result_pre{1}.w; result_pre{2}.w];
result_DD2 = matRad_fluenceOptimizationJO(dij_DD2,cst_DD,plnJO_DD2);



%% set pln
% only for protons, then only for photons and combine the plan in the end
% meta information for treatment plan (1) 
pln_DD2_robust(1).numOfFractions  = 5;
pln_DD2_robust(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_DD2_robust(1).machine         = 'Generic';

% beam geometry settings
pln_DD2_robust(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_DD2_robust(1).propStf.gantryAngles    = [-45 45];%[45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_DD2_robust(1).propStf.couchAngles     = zeros(numel(pln_DD2_robust(1).propStf.gantryAngles),1); % [?] ; 
pln_DD2_robust(1).propStf.numOfBeams      = numel(pln_DD2_robust(1).propStf.gantryAngles);
pln_DD2_robust(1).propStf.isoCenter       = ones(pln_DD2_robust(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_DD,ct,0);
% optimization settings
pln_DD2_robust(1).propDoseCalc.calcLET = 1;

pln_DD2_robust(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_DD2_robust(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_DD2_robust(1).propOpt.spatioTemp      = 0;
pln_DD2_robust(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_DD2_robust(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_DD2_robust(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_DD2_robust(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

% worst case scenario for robustness analysis
scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve scenarios for dose calculation and optimziation
pln_DD2_robust(1).multScen = matRad_multScen(ct,scenGenType);

% retrieve bio model parameters
pln_DD2_robust(1).bioParam = matRad_bioModel(pln_DD2_robust(1).radiationMode,quantityOpt, modelName);

plnJO_DD2_robust = matRad_plnWrapper(pln_DD2_robust);

%% multiple scenarios
modality = 'protons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_DD2_new';
multScen_DD2 = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_DD2(1:2),pln_DD2_robust,cst_DD,modality,[], saveDirectory);


%% LMO objectives
%cst{2,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,30));
cst_LxD = cst;
cst_LxD{1,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,0));
cst_LxD{4,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,0));
%% pln
% % meta information for treatment plan (1) 
pln_LxD1(1).numOfFractions  = 5;
pln_LxD1(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_LxD1(1).machine         = 'Generic';

% beam geometry settings
pln_LxD1(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_LxD1(1).propStf.gantryAngles    = 0;%[265 285]; %[90 180 270];% [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_LxD1(1).propStf.couchAngles     = zeros(numel(pln_LxD1(1).propStf.gantryAngles),1); % [?] ; 
pln_LxD1(1).propStf.numOfBeams      = numel(pln_LxD1(1).propStf.gantryAngles);
pln_LxD1(1).propStf.isoCenter       = ones(pln_LxD1(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_LxD,ct,0);
% optimization settings
pln_LxD1(1).propDoseCalc.calcLET = 1;

pln_LxD1(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_LxD1(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_LxD1(1).propOpt.spatioTemp      = 0;
pln_LxD1(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_LxD1(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_LxD1(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_LxD1(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
%quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
%=======================================> Model check error in bioModel
%modelName    = 'MCN';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions
modelName = 'constRBE';


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_LxD1(1).bioParam = matRad_bioModel(pln_LxD1(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_LxD1(1).multScen = matRad_multScen(ct,scenGenType);

%pln = pln(1);

% meta information for treatment plan (2) 
pln_LxD1(2).numOfFractions  = 25;
pln_LxD1(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_LxD1(2).machine         = 'Generic';

% beam geometry settings
pln_LxD1(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_LxD1(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_LxD1(2).propStf.couchAngles     = zeros(numel(pln_LxD1(2).propStf.gantryAngles),1);  % [?] ; 
pln_LxD1(2).propStf.numOfBeams      = numel(pln_LxD1(2).propStf.gantryAngles);
pln_LxD1(2).propStf.isoCenter       = ones(pln_LxD1(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_LxD,ct,0);
% optimization settings
pln_LxD1(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_LxD1(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_LxD1(2).propOpt.spatioTemp      = 0;
pln_LxD1(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_LxD1(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_LxD1(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_LxD1(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

%quantityOpt  = ['effect'];     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_LxD1(2).bioParam = matRad_bioModel(pln_LxD1(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_LxD1(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst_LxD = matRad_prepCst(cst_LxD, sparecst);

%pln = pln(2);

% Plan Wrapper
plnJO_LxD1 = matRad_plnWrapper(pln_LxD1);

%% Stf Wrapper
stf_LxD1 = matRad_stfWrapper(ct,cst_LxD,plnJO_LxD1);

%% Dij Calculation
dij_LxD1 = matRad_calcCombiDose(ct,stf_LxD1,plnJO_LxD1,cst_LxD,false);
% Dirty Dose Calculation
dij_LxD1 = matRad_calcDirtyDose(2,dij_LxD1,pln_LxD1);
dij_LxD1 = matRad_calcmLETDose(dij_LxD1,pln_LxD1);

%% compress optimized results into one cube 
%
dij_LxD1.precon = 1;
% dij.wInit = [result_pre{1}.w; result_pre{2}.w];
result_LxD1 = matRad_fluenceOptimizationJO(dij_LxD1,cst_LxD,plnJO_LxD1);



%% set pln
% only for protons, then only for photons and combine the plan in the end
% meta information for treatment plan (1) 
pln_LxD1_robust(1).numOfFractions  = 5;
pln_LxD1_robust(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_LxD1_robust(1).machine         = 'Generic';

% beam geometry settings
pln_LxD1_robust(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_LxD1_robust(1).propStf.gantryAngles    = 0;%[-45 45];%[45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_LxD1_robust(1).propStf.couchAngles     = zeros(numel(pln_LxD1_robust(1).propStf.gantryAngles),1); % [?] ; 
pln_LxD1_robust(1).propStf.numOfBeams      = numel(pln_LxD1_robust(1).propStf.gantryAngles);
pln_LxD1_robust(1).propStf.isoCenter       = ones(pln_LxD1_robust(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_LxD,ct,0);
% optimization settings
pln_LxD1_robust(1).propDoseCalc.calcLET = 1;

pln_LxD1_robust(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_LxD1_robust(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_LxD1_robust(1).propOpt.spatioTemp      = 0;
pln_LxD1_robust(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_LxD1_robust(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_LxD1_robust(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_LxD1_robust(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

% worst case scenario for robustness analysis
scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve scenarios for dose calculation and optimziation
pln_LxD1_robust(1).multScen = matRad_multScen(ct,scenGenType);

% retrieve bio model parameters
pln_LxD1_robust(1).bioParam = matRad_bioModel(pln_LxD1_robust(1).radiationMode,quantityOpt, modelName);

plnJO_LxD1_robust = matRad_plnWrapper(pln_LxD1_robust);

%% multiple scenarios
modality = 'protons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_LxD1_new';
multScen_LxD1 = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_LxD1(1),pln_LxD1_robust,cst_LxD,modality,[], saveDirectory);

%% pln
% % meta information for treatment plan (1) 
pln_LxD2(1).numOfFractions  = 5;
pln_LxD2(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_LxD2(1).machine         = 'Generic';

% beam geometry settings
pln_LxD2(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_LxD2(1).propStf.gantryAngles    = [-45 45];%[265 285]; %[90 180 270];% [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_LxD2(1).propStf.couchAngles     = zeros(numel(pln_LxD2(1).propStf.gantryAngles),1); % [?] ; 
pln_LxD2(1).propStf.numOfBeams      = numel(pln_LxD2(1).propStf.gantryAngles);
pln_LxD2(1).propStf.isoCenter       = ones(pln_LxD2(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_LxD,ct,0);
% optimization settings
pln_LxD2(1).propDoseCalc.calcLET = 1;

pln_LxD2(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_LxD2(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_LxD2(1).propOpt.spatioTemp      = 0;
pln_LxD2(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_LxD2(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_LxD2(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_LxD2(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
%quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
%=======================================> Model check error in bioModel
%modelName    = 'MCN';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions
modelName = 'constRBE';


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_LxD2(1).bioParam = matRad_bioModel(pln_LxD2(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_LxD2(1).multScen = matRad_multScen(ct,scenGenType);

%pln = pln(1);

% meta information for treatment plan (2) 
pln_LxD2(2).numOfFractions  = 25;
pln_LxD2(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_LxD2(2).machine         = 'Generic';

% beam geometry settings
pln_LxD2(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_LxD2(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_LxD2(2).propStf.couchAngles     = zeros(numel(pln_LxD2(2).propStf.gantryAngles),1);  % [?] ; 
pln_LxD2(2).propStf.numOfBeams      = numel(pln_LxD2(2).propStf.gantryAngles);
pln_LxD2(2).propStf.isoCenter       = ones(pln_LxD2(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_LxD,ct,0);
% optimization settings
pln_LxD2(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_LxD2(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_LxD2(2).propOpt.spatioTemp      = 0;
pln_LxD2(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_LxD2(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_LxD2(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_LxD2(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

%quantityOpt  = ['effect'];     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_LxD2(2).bioParam = matRad_bioModel(pln_LxD2(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_LxD2(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst_LxD = matRad_prepCst(cst_LxD, sparecst);

%pln = pln(2);

% Plan Wrapper
plnJO_LxD2 = matRad_plnWrapper(pln_LxD2);

%% Stf Wrapper
stf_LxD2 = matRad_stfWrapper(ct,cst_LxD,plnJO_LxD2);

%% Dij Calculation
dij_LxD2 = matRad_calcCombiDose(ct,stf_LxD2,plnJO_LxD2,cst_LxD,false);
% Dirty Dose Calculation
dij_LxD2 = matRad_calcDirtyDose(2,dij_LxD2,pln_LxD2);
dij_LxD2 = matRad_calcmLETDose(dij_LxD2,pln_LxD2);

%% compress optimized results into one cube 
%
dij_LxD2.precon = 1;
% dij.wInit = [result_pre{1}.w; result_pre{2}.w];
result_LxD2 = matRad_fluenceOptimizationJO(dij_LxD2,cst_LxD,plnJO_LxD2);



%% set pln
% only for protons, then only for photons and combine the plan in the end
% meta information for treatment plan (1) 
pln_LxD2_robust(1).numOfFractions  = 5;
pln_LxD2_robust(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_LxD2_robust(1).machine         = 'Generic';

% beam geometry settings
pln_LxD2_robust(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_LxD2_robust(1).propStf.gantryAngles    = [-45 45];%[45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_LxD2_robust(1).propStf.couchAngles     = zeros(numel(pln_LxD2_robust(1).propStf.gantryAngles),1); % [?] ; 
pln_LxD2_robust(1).propStf.numOfBeams      = numel(pln_LxD2_robust(1).propStf.gantryAngles);
pln_LxD2_robust(1).propStf.isoCenter       = ones(pln_LxD2_robust(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_LxD,ct,0);
% optimization settings
pln_LxD2_robust(1).propDoseCalc.calcLET = 1;

pln_LxD2_robust(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_LxD2_robust(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_LxD2_robust(1).propOpt.spatioTemp      = 0;
pln_LxD2_robust(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_LxD2_robust(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_LxD2_robust(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_LxD2_robust(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

% worst case scenario for robustness analysis
scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve scenarios for dose calculation and optimziation
pln_LxD2_robust(1).multScen = matRad_multScen(ct,scenGenType);

% retrieve bio model parameters
pln_LxD2_robust(1).bioParam = matRad_bioModel(pln_LxD2_robust(1).radiationMode,quantityOpt, modelName);

plnJO_LxD2_robust = matRad_plnWrapper(pln_LxD2_robust);

%% multiple scenarios
modality = 'protons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_LxD2_new';
multScen_LxD2 = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_LxD2(1:2),pln_LxD2_robust,cst_LxD,modality,[], saveDirectory);


%%
% meta information for treatment plan (2) 
pln_photonDD_robust.numOfFractions  = 25;
pln_photonDD_robust.radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_photonDD_robust.machine         = 'Generic';

% beam geometry settings
pln_photonDD_robust.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_photonDD_robust.propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_photonDD_robust.propStf.couchAngles     = zeros(numel(pln_photonDD_robust.propStf.gantryAngles),1);  % [?] ; 
pln_photonDD_robust.propStf.numOfBeams      = numel(pln_photonDD_robust.propStf.gantryAngles);
pln_photonDD_robust.propStf.isoCenter       = ones(pln_photonDD_robust.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_DD,ct,0);
% optimization settings
pln_photonDD_robust.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_photonDD_robust.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_photonDD_robust.propOpt.spatioTemp      = 0;
pln_photonDD_robust.propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_photonDD_robust.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_photonDD_robust.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_photonDD_robust.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_photonDD_robust.bioParam = matRad_bioModel(pln_photonDD_robust.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_photonDD_robust.multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  

%% multiple scenarios
modality = 'photons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Photons_DD_new';
multScen_DDphotons = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_DD1(2:end),pln_photonDD_robust,cst_DD,modality,[], saveDirectory);

%%
% meta information for treatment plan (2) 
pln_photonLxD_robust.numOfFractions  = 25;
pln_photonLxD_robust.radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_photonLxD_robust.machine         = 'Generic';

% beam geometry settings
pln_photonLxD_robust.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_photonLxD_robust.propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_photonLxD_robust.propStf.couchAngles     = zeros(numel(pln_photonLxD_robust.propStf.gantryAngles),1);  % [?] ; 
pln_photonLxD_robust.propStf.numOfBeams      = numel(pln_photonLxD_robust.propStf.gantryAngles);
pln_photonLxD_robust.propStf.isoCenter       = ones(pln_photonLxD_robust.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_LxD,ct,0);
% optimization settings
pln_photonLxD_robust.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_photonLxD_robust.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_photonLxD_robust.propOpt.spatioTemp      = 0;
pln_photonLxD_robust.propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_photonLxD_robust.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_photonLxD_robust.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_photonLxD_robust.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_photonLxD_robust.bioParam = matRad_bioModel(pln_photonLxD_robust.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_photonLxD_robust.multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  

%% multiple scenarios
modality = 'photons';
saveDirectory = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Photons_LxD_new';
multScen_LxDphotons = matRad_calcParticleDoseMultipleScenariosExtended(ct,stf_LxD1(2:end),pln_photonLxD_robust,cst_LxD,modality,[], saveDirectory);

%% get meta information
% Photons 9 beams DD
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Photons_DD_new';
scenIdxToInclude = 1:9;
meta_photonsDD = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Photons 9 beams LxD
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Photons_LxD_new';
scenIdxToInclude = 1:9;
meta_photonsLxD = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons DD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_DD1_new';
scenIdxToInclude = 1:9;
meta_DD1 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons DD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_DD2_new';
scenIdxToInclude = 1:9;
meta_DD2 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons LxD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_LxD1_new';
scenIdxToInclude = 1:9;
meta_LxD1 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);

% Protons LxD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_LxD2_new';
scenIdxToInclude = 1:9;
meta_LxD2 = matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude);


%% load meta information

% photons 9 beams DD
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Photons_DD_new';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Photons_DD_new\dijTemplate.mat")
dij_photonsDD = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    Scen_photonsDD{scen} = matRad_loadScenariosFromMeta(saveDir, meta_photonsDD(scen), dij_photonsDD, 1); 
end

% photons 9 beams LxD
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Photons_LxD_new';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Photons_LxD_new\dijTemplate.mat")
dij_photonsLxD = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    Scen_photonsLxD{scen} = matRad_loadScenariosFromMeta(saveDir, meta_photonsLxD(scen), dij_photonsLxD, 1); 
end


% protons DD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_DD1_new';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_DD1_new\dijTemplate.mat")
dij_protonsDD1 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    Scen_protonsDD1{scen} = matRad_loadScenariosFromMeta(saveDir, meta_DD1(scen), dij_protonsDD1, 1); 
end

% protons DD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_DD2_new';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_DD2_new\dijTemplate.mat")
dij_protonsDD2 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    Scen_protonsDD2{scen} = matRad_loadScenariosFromMeta(saveDir, meta_DD2(scen), dij_protonsDD2, 1); 
end

% protons LxD 0°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_LxD1_new';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_LxD1_new\dijTemplate.mat")
dij_protonsLxD1 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    Scen_protonsLxD1{scen} = matRad_loadScenariosFromMeta(saveDir, meta_LxD1(scen), dij_protonsLxD1, 1); 
end

% protons LxD -45° and 45°
saveDir = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_LxD2_new';
load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\RobustScenDD_new\Protons_LxD2_new\dijTemplate.mat")
dij_protonsLxD2 = dijTemplate;
stringLength = 0;
nScenarios = size(scenIdxToInclude,2);
for scen = scenIdxToInclude
    fprintf(repmat('\b',1,stringLength));
    stringLength = fprintf('\tLoading scenario: %u/%u\n',scen,nScenarios);
    Scen_protonsLxD2{scen} = matRad_loadScenariosFromMeta(saveDir, meta_LxD2(scen), dij_protonsLxD2, 1); 
end

%% Weighting physicalDose

for scen = scenIdxToInclude
    photonDDPhysD{scen} = Scen_photonsDD{1,scen}.physicalDose{1,1} * result_DD1{1,2}.w;
    photonLxDPhysD{scen} = Scen_photonsLxD{1,scen}.physicalDose{1,1} * result_LxD1{1,2}.w;
    protonDD1PhysD{scen} = Scen_protonsDD1{1,scen}.physicalDose{1,1} * result_DD1{1,1}.w;
    protonDD2PhysD{scen} = Scen_protonsDD2{1,scen}.physicalDose{1,1} * result_DD2{1,1}.w;
    protonLxD1PhysD{scen} = Scen_protonsLxD1{1,scen}.physicalDose{1,1} * result_LxD1{1,1}.w;
    protonLxD2PhysD{scen} = Scen_protonsLxD2{1,scen}.physicalDose{1,1} * result_LxD2{1,1}.w;
end

%% Reshaping

for scen = scenIdxToInclude
    PhotonDDPhysD{scen} = reshape(photonDDPhysD{scen}, 167, 167, 129);
    PhotonLxDPhysD{scen} = reshape(photonLxDPhysD{scen}, 167, 167, 129);
    ProtonDD1PhysD{scen} = reshape(protonDD1PhysD{scen}, 167, 167, 129);
    ProtonDD2PhysD{scen} = reshape(protonDD2PhysD{scen}, 167, 167, 129);
    ProtonLxD1PhysD{scen} = reshape(protonLxD1PhysD{scen}, 167, 167, 129);
    ProtonLxD2PhysD{scen} = reshape(protonLxD2PhysD{scen}, 167, 167, 129);

end

%% Combining
for scen = scenIdxToInclude
    CombiDD1_sameScen{scen} = 5 * ProtonDD1PhysD{scen} + 25 * PhotonDDPhysD{scen};
    CombiDD2_sameScen{scen} = 5 * ProtonDD2PhysD{scen} + 25 * PhotonDDPhysD{scen};
    CombiLxD1_sameScen{scen} = 5 * ProtonLxD1PhysD{scen} + 25 * PhotonLxDPhysD{scen};
    CombiLxD2_sameScen{scen} = 5 * ProtonLxD2PhysD{scen} + 25 * PhotonLxDPhysD{scen};
end

%% DVH
physDD = 5 * result_DD1{1,1}.physicalDose + 25 * result_DD1{1,2}.physicalDose;
physDD2 = 5 * result_DD2{1,1}.physicalDose + 25 * result_DD2{1,2}.physicalDose;

physLxD = 5 * result_LxD1{1,1}.physicalDose + 25 * result_LxD1{1,2}.physicalDose;
physLxD2 = 5 * result_LxD2{1,1}.physicalDose + 25 * result_LxD2{1,2}.physicalDose;
 
DVH_DD1_opt = matRad_calcDVH(cst_DD,physDD);
DVH_DD2_opt = matRad_calcDVH(cst_DD,physDD2);

DVH_LxD1_opt = matRad_calcDVH(cst_LxD,physLxD);
DVH_LxD2_opt = matRad_calcDVH(cst_LxD,physLxD2);

%% plot

%% Plotting DVH
color = ['b', 'm', 'g' , 'r'];

figure
structure = [1 2 3];
for i = [1 2 3]
    plot(DVH_DD1_opt(structure(i)).doseGrid,DVH_DD1_opt(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',2)
    hold on
    plot(DVH_DD2_opt(structure(i)).doseGrid,DVH_DD2_opt(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',2)
    hold on
    plot(DVH_LxD1_opt(structure(i)).doseGrid,DVH_LxD1_opt(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',2)
    hold on
    plot(DVH_LxD2_opt(structure(i)).doseGrid,DVH_LxD2_opt(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',2)
    hold on    
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
    plot(DVH_DD1{4}(structure(i)).doseGrid,DVH_DD1{4}(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',2)
    hold on
    plot(DVH_DD2{4}(structure(i)).doseGrid,DVH_DD2{4}(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',2)
    hold on
    plot(DVH_LxD1{4}(structure(i)).doseGrid,DVH_LxD1{4}(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',2)
    hold on
    plot(DVH_LxD2{4}(structure(i)).doseGrid,DVH_LxD2{4}(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',2)
    hold on    
end
