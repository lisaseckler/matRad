matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-06;
load 'NEW-BOXPHANTOM-Overlap.mat'

% changing alphaX
cst{2,5}.alphaX = 0.5;
% cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(1000,30));
% cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(20,0));
% cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));
% cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));
% cst{5,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(20,0));
cst{1,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
%cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,30));
%cst{1,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,120));

%% add margin
% cube = zeros(ct.cubeDim);
% cube(cst{1,4}{1}) = 1;
% vResolution = ct.resolution;
% vMargin = [];
% vMargin.x = 5;
% vMargin.y = 5;
% vMargin.z = 5;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin,1);
% 
% cst{4,1}    = 3;
% cst{4,2}    = 'Margin';
% cst{4,3}    = 'OAR';
% cst{4,4}{1} = find(mVOIEnlarged);
% cst{4,5}    = cst{1,5};
% 
% cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40)); 

%%
% meta information for treatment plan (1) 
pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
% pln(1).propStf.gantryAngles    = [ 45 0 -45]; % [?] ;
pln(1).propStf.gantryAngles    = [90];
pln(1).propStf.couchAngles     = zeros(numel(pln(1).propStf.gantryAngles),1); % [?] ; 
pln(1).propStf.numOfBeams      = numel(pln(1).propStf.gantryAngles);
pln(1).propStf.isoCenter       = ones(pln(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(1).propDoseCalc.calcLET = 1;

pln(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).propOpt.spatioTemp      = 0;
pln(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(1).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(1).bioParam = matRad_bioModel(pln(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(1).multScen = matRad_multScen(ct,scenGenType);
% 
% meta information for treatment plan (2) 
pln(2).numOfFractions  = 25;
pln(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [0:72:359]; % [?] ;
pln(2).propStf.couchAngles     = zeros(numel(pln(2).propStf.gantryAngles),1);  % [?] ; 
pln(2).propStf.numOfBeams      = numel(pln(2).propStf.gantryAngles);
pln(2).propStf.isoCenter       = ones(pln(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(2).propOpt.spatioTemp      = 0;
pln(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(2).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(2).bioParam = matRad_bioModel(pln(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);
% Plan Wrapper
plnJO = matRad_plnWrapper(pln);
% Stf Wrapper
stf = matRad_stfWrapper(ct,cst,plnJO);
% Dij Calculation
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
% Dirty Dose Calculation
dij = matRad_calcDirtyDose(2,dij,pln);
% mLETDose sumed up
%dij = matRad_calcmLETDose(dij,pln);
% Fluence optimization 
result = matRad_fluenceOptimizationJO(dij,cst,plnJO);
resultGUI = matRad_calcResultGUIstruct(result);
% pln_original = pln;
% pln = pln(2);
% matRadGUI

%% DVH calculation


%% Visualization
% slice = 59;
% 
% photon_plan = resultGUI{2};
% proton_plan = resultGUI{1};
% totalPlan = pln(1).numOfFractions.*proton_plan.(quantityOpt) + pln(2).numOfFractions.*photon_plan.(quantityOpt);
% 
% f = figure;
% subplot(1,3,1);
%     imagesc(proton_plan.(quantityOpt)(:,:,slice));
%     matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
%     title('Proton Plan');
% subplot(1,3,2);
%     imagesc(photon_plan.(quantityOpt)(:,:,slice));
%     matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
%     title('Photon Plan');
% subplot(1,3,3);
%     imagesc(totalPlan(:,:,slice));
%     matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
%     title('Total Plan');