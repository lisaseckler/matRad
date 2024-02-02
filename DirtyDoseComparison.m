% matRad_rc;
% matRad_cfg = MatRad_Config.instance();
% matRad_cfg.propOpt.defaultMaxIter = 50000;
% matRad_cfg.propOpt.defaultAccChangeTol = 1e-04;
load NEW-BOXPHANTOM-Overlap.mat

% changing alphaX
%cst{2,5}.alphaX = 0.1;

cst{1,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(800,60));
%cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));
%cst{3,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosing(100,0));

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
result = matRad_fluenceOptimizationJO(dij,cst,plnJO);
resultGUI_1 = matRad_calcResultGUIstruct(result);
%result_first = result;
%dij_first = dij;
result_without = result;

save("MixedModalities_withoutDD.mat","-v7.3")
%% with DD in OAR2
clear
load 'NEW-BOXPHANTOM-Overlap.mat'

% changing alphaX
%cst{2,5}.alphaX = 0.1;

cst{1,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(800,60));
cst{3,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,30));

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
result = matRad_fluenceOptimizationJO(dij,cst,plnJO);
resultGUI = matRad_calcResultGUIstruct(result);
% cst_secondD = cst;
% dij_secondD = dij2D;
% result_secondD = result2D;
% resultGUI_secondD = result2D{1,1};
result_DD = result;

save("MixedModalities_withDD.mat","-v7.3")

%%
clear
load 'NEW-BOXPHANTOM-Overlap.mat'

% changing alphaX
%cst{2,5}.alphaX = 0.5;

%cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));
cst{1,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(800,60));
cst{3,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,120));
%cst_third = cst;

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
result = matRad_fluenceOptimizationJO(dij,cst,plnJO);
resultGUI = matRad_calcResultGUIstruct(result);
% dij_third = dij;
% resultGUI_third = result3{1,1};
result_mL = result;
save("MixedModalities_withmL.mat","-v7.3")
%% New Case
% clear("cst","ct");
% load 'NEW-BOXPHANTOM-Overlap.mat'
% 
% % changing alphaX
% cst{2,5}.alphaX = 0.5;
% 
% cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));
% cst{3,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,30));
% 
% cst{1,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
% cst_fourthD = cst;
% 
% % meta information for treatment plan (1) 
% pln4D(1).numOfFractions  = 5;
% pln4D(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
% pln4D(1).machine         = 'Generic';
% 
% % beam geometry settings
% pln4D(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
% % pln(1).propStf.gantryAngles    = [ 45 0 -45]; % [?] ;
% pln4D(1).propStf.gantryAngles    = [90];
% pln4D(1).propStf.couchAngles     = zeros(numel(pln4D(1).propStf.gantryAngles),1); % [?] ; 
% pln4D(1).propStf.numOfBeams      = numel(pln4D(1).propStf.gantryAngles);
% pln4D(1).propStf.isoCenter       = ones(pln4D(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% % optimization settings
% pln4D(1).propDoseCalc.calcLET = 1;
% 
% pln4D(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
% pln4D(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
% pln4D(1).propOpt.spatioTemp      = 0;
% pln4D(1).propOpt.STscenarios     = 2;
% %pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)
% 
% % dose calculation settings
% pln4D(1).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
% pln4D(1).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
% pln4D(1).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% % pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
% quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
% %=======================================> Model check error in bioModel
% modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
%                                    % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
%                                    % LEM: Local Effect Model for carbon ions
% 
% 
% scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          
% 
% % retrieve bio model parameters
% pln4D(1).bioParam = matRad_bioModel(pln4D(1).radiationMode,quantityOpt, modelName);
% 
% % retrieve scenarios for dose calculation and optimziation
% pln4D(1).multScen = matRad_multScen(ct,scenGenType);
% % 
% % meta information for treatment plan (2) 
% pln4D(2).numOfFractions  = 25;
% pln4D(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
% pln4D(2).machine         = 'Generic';
% 
% % beam geometry settings
% pln4D(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
% pln4D(2).propStf.gantryAngles    = [0:72:359]; % [?] ;
% pln4D(2).propStf.couchAngles     = zeros(numel(pln4D(2).propStf.gantryAngles),1);  % [?] ; 
% pln4D(2).propStf.numOfBeams      = numel(pln4D(2).propStf.gantryAngles);
% pln4D(2).propStf.isoCenter       = ones(pln4D(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% % optimization settings
% pln4D(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
% pln4D(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
% pln4D(2).propOpt.spatioTemp      = 0;
% pln4D(2).propOpt.STscenarios     = 5;
% %pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)
% 
% % dose calculation settings
% pln4D(2).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
% pln4D(2).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
% pln4D(2).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% % pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;
% 
% quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
% modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
%                                    % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
%                                    % LEM: Local Effect Model for carbon ions
% 
% 
% scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          
% 
% % retrieve bio model parameters
% pln4D(2).bioParam = matRad_bioModel(pln4D(2).radiationMode,quantityOpt, modelName);
% 
% % retrieve scenarios for dose calculation and optimziation
% pln4D(2).multScen = matRad_multScen(ct,scenGenType);
% 
% % prepping cst 
% % placing alpha/beta ratios in cst{:,6},
% % different alpha beta ration for each obj of a structure  
% sparecst = 0;
% 
% cst = matRad_prepCst(cst, sparecst);
% % Plan Wrapper
% plnJO4D = matRad_plnWrapper(pln4D);
% % Stf Wrapper
% stf4D = matRad_stfWrapper(ct,cst,plnJO4D);
% % Dij Calculation
% dij4D = matRad_calcCombiDose(ct,stf4D,plnJO4D,cst,false);
% % Dirty Dose Calculation
% dij4D = matRad_calcDirtyDose(2,dij4D,pln4D);
% result4D = matRad_fluenceOptimizationJO(dij4D,cst,plnJO4D);
% resultGUI4D = matRad_calcResultGUIstruct(result4D);
% resultGUI_fourthD = result4D{1,1};
% 
% 
% %% DVH Calculation
% 
% dvh_withoutDD = matRad_calcDVH(cst_first,resultGUI_first.dirtyDose);
% dvh_withoutDD_mL = matRad_calcDVH(cst_first,resultGUI_first.mLETDose);
% dvh_withDD = matRad_calcDVH(cst_second,resultGUI_second.dirtyDose);
% dvh_withDD_mL = matRad_calcDVH(cst_second,resultGUI_second.mLETDose);
% 
% dvh_withDD_D = matRad_calcDVH(cst_secondD,resultGUI_secondD.dirtyDose);
% dvh_withDD_D_05 = matRad_calcDVH(cst_fourthD,resultGUI_fourthD.dirtyDose);
% dvh_withDD_D_mL = matRad_calcDVH(cst_secondD,resultGUI_secondD.mLETDose);
% dvh_withDD_D_mL_05 = matRad_calcDVH(cst_fourthD,resultGUI_fourthD.mLETDose);
% 
% dvh_withoutDD_05 = matRad_calcDVH(cst_third,resultGUI_third.dirtyDose);
% dvh_withoutDD_mL_05 = matRad_calcDVH(cst_third,resultGUI_third.mLETDose);
% dvh_withDD_05 = matRad_calcDVH(cst_fourth,resultGUI_fourth.dirtyDose);
% dvh_withDD_mL_05 = matRad_calcDVH(cst_fourth,resultGUI_fourth.mLETDose);
% 
% figure
% subplot(2,1,1)
% for i = [1 2 3]
%     plot(dvh_withoutDD(i).doseGrid,dvh_withoutDD(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
%     hold on
%     plot(dvh_withoutDD_05(i).doseGrid,dvh_withoutDD_05(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
%     hold on
%     plot(dvh_withDD(i).doseGrid,dvh_withDD(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-.','LineWidth',1.5)
%     hold on
%     plot(dvh_withDD_05(i).doseGrid,dvh_withDD_05(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle',':','LineWidth',1.5)
% end
% title('alphaX 0.1')
% xlabel('dirty dose in Gy')
% ylabel('Volume in %')
% legend('withoutDD_01','withDD_01','withoutDD_05','withDD_05')
% subplot(2,1,2)
% for i = [1 2 3]
%     plot(dvh_withoutDD_mL_05(i).doseGrid,dvh_withoutDD_mL_05(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
%     hold on
%     plot(dvh_withoutDD_mL(i).doseGrid,dvh_withoutDD_mL(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
%     hold on
%     plot(dvh_withDD_mL_05(i).doseGrid,dvh_withDD_mL_05(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-.','LineWidth',1.5)
%     hold on
%     plot(dvh_withDD_mL(i).doseGrid,dvh_withDD_mL(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle',':','LineWidth',1.5)
% end
% title('alphaX 0.1')
% xlabel('dirty dose in Gy')
% ylabel('Volume in %')
% legend('withoutDD_01','withDD_01','withoutDD_05','withDD_05')
% 
% figure
% subplot(2,1,1)
% for i = [1 2 3]
%     plot(dvh_withDD_D(i).doseGrid,dvh_withDD_D(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
%     hold on
%     plot(dvh_withDD_D_05(i).doseGrid,dvh_withDD_D_05(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
% end
% xlabel('dirty dose in Gy')
% ylabel('Volume in %')
% legend('withoutDD','withDD_D')
% title('DVH: PTV violett, OAR red, Body brown')
% subplot(2,1,2)
% for i = [1 2 3]
%     plot(dvh_withDD_D_mL(i).doseGrid,dvh_withDD_D_mL(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
%     hold on
%     plot(dvh_withDD_D_mL_05(i).doseGrid,dvh_withDD_D_mL_05(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
% end
% xlabel('mLETDose')
% ylabel('Volume in %')
% legend('withDD_mL','withDD_D_mL')
% 
% figure
% subplot(2,1,1)
% for i = [1 2 3]
%     plot(dvh_withDD(i).doseGrid,dvh_withDD(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
%     hold on
%     plot(dvh_withDD_D(i).doseGrid,dvh_withDD_D(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
% end
% xlabel('dirty dose in Gy')
% ylabel('Volume in %')
% legend('withoutDD','withDD_D')
% title('DVH: PTV violett, OAR red, Body brown')
% subplot(2,1,2)
% for i = [1 2 3]
%     plot(dvh_withDD_mL(i).doseGrid,dvh_withDD_mL(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
%     hold on
%     plot(dvh_withDD_D_mL(i).doseGrid,dvh_withDD_D_mL(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
% end
% xlabel('mLETDose')
% ylabel('Volume in %')
% legend('withDD_mL','withDD_D_mL')