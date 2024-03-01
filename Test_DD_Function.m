matRad_rc
clear
load("NEW-BOXPHANTOM-Overlap.mat")

cst{3,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,120));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(800,60));
cst_N_ml = cst;

plnN_ml(1).numOfFractions  = 30;
plnN_ml(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
plnN_ml(1).machine         = 'Generic';

% beam geometry settings
plnN_ml(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
% pln(1).propStf.gantryAngles    = [ 45 0 -45]; % [?] ;
plnN_ml(1).propStf.gantryAngles    = [90];
plnN_ml(1).propStf.couchAngles     = zeros(numel(plnN_ml(1).propStf.gantryAngles),1); % [?] ; 
plnN_ml(1).propStf.numOfBeams      = numel(plnN_ml(1).propStf.gantryAngles);
plnN_ml(1).propStf.isoCenter       = ones(plnN_ml(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
plnN_ml(1).propDoseCalc.calcLET = 1;

plnN_ml(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
plnN_ml(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
plnN_ml(1).propOpt.spatioTemp      = 0;
plnN_ml(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
plnN_ml(1).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
plnN_ml(1).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
plnN_ml(1).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
plnN_ml(1).bioParam = matRad_bioModel(plnN_ml(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
plnN_ml(1).multScen = matRad_multScen(ct,scenGenType);

sparecst = 0;

cst_N_ml = matRad_prepCst(cst_N_ml, sparecst);

plnJON_ml = matRad_plnWrapper(plnN_ml);

stfN_ml = matRad_stfWrapper(ct,cst_N_ml,plnJON_ml);
% Dij Calculation
dijN_ml = matRad_calcCombiDose(ct,stfN_ml,plnJON_ml,cst_N_ml,false);
%dijN = matRad_calcParticleDose(ct,stfN,plnN,cst);
% Dirty Dose Calculation
dijN_ml = matRad_calcDirtyDose(2,dijN_ml,plnN_ml);
%%
resultN_ml = matRad_fluenceOptimizationJO(dijN_ml,cst_N_ml,plnJON_ml);
%%
%resultN = matRad_fluenceOptimization(dijN,cst,plnN);

cst_N_DD = cst;
dij_N_DD = dijN;
pln_N_DD = plnN;
result_N_DD = resultN;
stf_N_DD = stfN;

save("mix_DD4.mat","cst_N_DD","dij_N_DD","pln_N_DD","result_N_DD","stf_N_DD","-v7.3")
%%
clear
load("NEW-BOXPHANTOM-Overlap.mat")

cst{3,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,120));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(800,60));

plnN(1).numOfFractions  = 30;
plnN(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
plnN(1).machine         = 'Generic';

% beam geometry settings
plnN(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
% pln(1).propStf.gantryAngles    = [ 45 0 -45]; % [?] ;
plnN(1).propStf.gantryAngles    = [90];
plnN(1).propStf.couchAngles     = zeros(numel(plnN(1).propStf.gantryAngles),1); % [?] ; 
plnN(1).propStf.numOfBeams      = numel(plnN(1).propStf.gantryAngles);
plnN(1).propStf.isoCenter       = ones(plnN(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
plnN(1).propDoseCalc.calcLET = 1;

plnN(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
plnN(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
plnN(1).propOpt.spatioTemp      = 0;
plnN(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
plnN(1).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
plnN(1).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
plnN(1).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
plnN(1).bioParam = matRad_bioModel(plnN(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
plnN(1).multScen = matRad_multScen(ct,scenGenType);

sparecst = 0;

cst = matRad_prepCst(cst, sparecst);

plnJON = matRad_plnWrapper(plnN);

stfN = matRad_stfWrapper(ct,cst,plnJON);
% Dij Calculation
dijN = matRad_calcCombiDose(ct,stfN,plnJON,cst,false);
%dijN = matRad_calcParticleDose(ct,stfN,plnN,cst);
% Dirty Dose Calculation
dijN = matRad_calcDirtyDose(2,dijN,plnN);
resultN = matRad_fluenceOptimizationJO(dijN,cst,plnJON);
%resultN = matRad_fluenceOptimization(dijN,cst,plnN);

cst_N_mL = cst;
dij_N_mL = dijN;
pln_N_mL = plnN;
result_N_mL = resultN;
stf_N_mL = stfN;

save("mix_mL4.mat","cst_N_mL","dij_N_mL","pln_N_mL","result_N_mL","stf_N_mL","-v7.3")

%% Visualisation

% cube = result_O_DD.RBExD;
% plane = 3;
% slice = 80;
% doseWindow = [0 max(cube(:))];
% 
% figure,
% subplot(2,2,1)
% matRad_plotSliceWrapper(gca,ct,cst_O_DD,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('Old branch DD, RBExD')
% zoom(1.3)
% 
% cube = result_O_mL.RBExD;
% plane = 3;
% slice = 80;
% doseWindow = [0 max(cube(:))];
% 
% subplot(2,2,2)
% matRad_plotSliceWrapper(gca,ct,cst_O_mL,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('Old branch mL, RBExD')
% zoom(1.3)
% 
% cube = resultN{1}.RBExD;
% plane = 3;
% slice = 80;
% doseWindow = [0 max(cube(:))];
% 
% subplot(2,2,3)
% matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('New branch DD, RBExD')
% zoom(1.3)
% % 
% cube = resultN_ml{1}.RBExD;
% plane = 3;
% slice = 80;
% doseWindow = [0 max(cube(:))];
% 
% subplot(2,2,4)
% matRad_plotSliceWrapper(gca,ct,cst_N_ml,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('New branch mL, RBExD')
% zoom(1.3)
% 
% cube = result_O_DD.RBExD - resultN{1}.RBExD;
% plane = 3;
% slice = 80;
% doseWindow = [min(cube(:)) max(cube(:))];
% 
% figure
% subplot(2,1,1)
% matRad_plotSliceWrapper(gca,ct,cst_O_DD,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('Difference DDs')
% zoom(1.3)
% 
% cube = resultN{1}.RBExD - resultN_ml{1}.RBExD;
% plane = 3;
% slice = 80;
% doseWindow = [min(cube(:)) max(cube(:))];
% 
% subplot(2,1,2)
% matRad_plotSliceWrapper(gca,ct,cst_O_DD,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('Difference DD and mL')
% zoom(1.3)