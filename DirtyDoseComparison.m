matRad_rc;
clear
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-06;
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

cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30)); 
%cst{4,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));


% changing alphaX
cst{2,5}.alphaX = 0.1;
%%
cst{1,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));

cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,60));
% cst{2,6}{3} = struct(DoseObjectives.matRad_SquaredUnderdosing(100,60));
%cst{2,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(100,61));
%cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));
%cst{3,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosing(100,0));

%%

% meta information for treatment plan (1) 
pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [90 180]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
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
pln(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(1).bioParam = matRad_bioModel(pln(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(1).multScen = matRad_multScen(ct,scenGenType);

%pln = pln(1);

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
pln(2).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
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
cst_without = cst;
% Plan Wrapper
plnJO = matRad_plnWrapper(pln);
% Stf Wrapper
stf = matRad_stfWrapper(ct,cst,plnJO);
% Dij Calculation
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
% Dirty Dose Calculation
dij = matRad_calcDirtyDose(2,dij,pln);
result_without = matRad_fluenceOptimizationJO(dij,cst,plnJO);
% resultGUI_1 = matRad_calcResultGUIstruct(result);
%result_first = result;
%dij_first = dij;
%result_without = result;

physDose_without = result_without{1,1}.physicalDose * 5 + result_without{1,2}.physicalDose * 25;
RBExD_without = result_without{1,1}.physicalDose * 1.1 * 5 + result_without{1,2}.physicalDose * 1.1 * 25;
%effect_without = result_without{1,1}.effect * 5 + result_without{1,2}.effect * 25;
dirtyDose_without = result_without{1,1}.dirtyDose * 5 + result_without{1,2}.dirtyDose * 25;
%LET_without = result_without{1,1}.LET * 5 + result_without{1,2}.LET * 5;
LET_without = (result_without{1,1}.LET .* result_without{1,1}.physicalDose * 5)./physDose_without;
%save("ProtonPhoton_TG119_MixedModalities_without.mat","-v7.3")

cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));
cst{4,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));
cst = matRad_prepCst(cst, sparecst);
cst_DD = cst;
stf = matRad_stfWrapper(ct,cst,plnJO);
% Dij Calculation
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
% Dirty Dose Calculation
dij = matRad_calcDirtyDose(2,dij,pln);
result_DD = matRad_fluenceOptimizationJO(dij,cst,plnJO);
% resultGUI_1 = matRad_calcResultGUIstruct(result);
%result_first = result;
%dij_first = dij;
%result_DD = result;
%
physDose_DD = result_DD{1,1}.physicalDose * 5 + result_DD{1,2}.physicalDose * 25;
RBExD_DD = result_DD{1,1}.physicalDose * 1.1 * 5 + result_DD{1,2}.physicalDose * 1.1 * 25;
%effect = result_DD{1,1}.effect * 5 + result_DD{1,2}.effect * 25;
dirtyDose_DD = result_DD{1,1}.dirtyDose * 5 + result_DD{1,2}.dirtyDose * 25;
%LET_DD = result_DD{1,1}.LET * 5 + result_DD{1,2}.LET * 5;
LET_DD = (result_DD{1,1}.LET .* result_DD{1,1}.physicalDose * 5)./physDose_DD;
%save("ProtonPhoton_TG119_MixedModalities_DD.mat","-v7.3")


cst{1,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,0));
cst{4,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,0));
cst = matRad_prepCst(cst, sparecst);
cst_mL = cst;
stf = matRad_stfWrapper(ct,cst,plnJO);
% Dij Calculation
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
% Dirty Dose Calculation
dij = matRad_calcDirtyDose(2,dij,pln);
dij = matRad_calcmLETDose(dij,pln);
result_mL2 = matRad_fluenceOptimizationJO(dij,cst,plnJO);
% resultGUI_1 = matRad_calcResultGUIstruct(result);
%result_first = result;
%dij_first = dij;
%result_mL = result;
%
physDose_mL2 = result_mL2{1,1}.physicalDose * 5 + result_mL2{1,2}.physicalDose * 25;
RBExD_mL2 = result_mL2{1,1}.physicalDose * 1.1 * 5 + result_mL2{1,2}.physicalDose * 1.1 * 25;
%effect = result_mL2{1,1}.effect * 5 + result_mL2{1,2}.effect * 25;
dirtyDose_mL2 = result_mL2{1,1}.dirtyDose * 5 + result_mL2{1,2}.dirtyDose * 25;
%LET_mL = result_mL{1,1}.LET * 5 + result_mL{1,2}.LET * 5;
LET_mL2 = (result_mL2{1,1}.LET .* result_mL2{1,1}.physicalDose * 5)./physDose_mL2;
%save("ProtonPhoton_TG119_MixedModalities_mL.mat","-v7.3")
%% BED
% BED = D*(1+d/(alpha/beta))
BED_without_proton = result_without{1,1}.physicalDose * pln(1).numOfFractions .*(1 + result_without{1,1}.physicalDose./(0.1./0.05));
BED_without_proton2 = result_without{1,1}.effect./0.1;
BED_DD_proton = result_DD{1,1}.physicalDose * pln(1).numOfFractions .*(1 + result_DD{1,1}.physicalDose./(0.1./0.05));
BED_DD_proton2 = result_DD{1,1}.effect./0.1;
BED_mL2_proton = result_mL2{1,1}.physicalDose * pln(1).numOfFractions .*(1 + result_mL2{1,1}.physicalDose./(0.1./0.05));
BED_mL_proton2 = result_mL{1,1}.effect./0.1;

BED_without_photon = result_without{1,2}.physicalDose * pln(2).numOfFractions .*(1 + result_without{1,2}.physicalDose./(0.1./0.05));
BED_without_photon2 = result_without{1,2}.effect./0.1;
BED_DD_photon = result_DD{1,2}.physicalDose * pln(2).numOfFractions .*(1 + result_DD{1,2}.physicalDose./(0.1./0.05));
BED_DD_photon2 = result_DD{1,2}.effect./0.1;
BED_mL2_photon = result_mL2{1,2}.physicalDose * pln(2).numOfFractions .*(1 + result_mL2{1,2}.physicalDose./(0.1./0.05));
BED_mL_photon2 = result_mL{1,2}.effect./0.1;

BED_without_combi = BED_without_proton + BED_without_photon;
BED_without_combi2 = BED_without_proton2 + BED_without_photon2;
BED_DD_combi = BED_DD_proton + BED_DD_photon;
BED_DD_combi2 = BED_DD_proton2 + BED_DD_photon2;
BED_mL2_combi = BED_mL2_proton + BED_mL2_photon;
BED_mL_combi2 = BED_mL_proton2 + BED_mL_photon2;


%% Visualisation
cube = effect;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total physicalDose with mL objective body 20')
zoom(1.3)

%% subplot total

cube = physDose_without;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total physicalDose without objective')
zoom(1.3)

cube = physDose_DD;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total physicalDose with DD objective')
zoom(1.3)

cube = physDose_mL2;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total physicalDose with mL objective')
zoom(1.3)

cube = dirtyDose_without;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total dirtyDose without objective')
zoom(1.3)

cube = dirtyDose_DD;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total dirtyDose with DD objective')
zoom(1.3)

cube = dirtyDose_mL2;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total dirtyDose with mL objective')
zoom(1.3)

cube = LET_without;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total LET without objective')
zoom(1.3)

cube = LET_DD;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total LET with DD objective')
zoom(1.3)

cube = LET_mL2;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total LET with mL objective')
zoom(1.3)

%% subplot only physicalDose

cube = result_without{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('proton physicalDose without objective')
zoom(1.3)

cube = result_without{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photon physicalDose without objective')
zoom(1.3)

cube = physDose_without;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total physicalDose without objective')
zoom(1.3)

cube = result_DD{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('proton physicalDose with DD objective')
zoom(1.3)

cube = result_DD{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photon physicalDose with DD objective')
zoom(1.3)

cube = physDose_DD;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total physicalDose with DD objective')
zoom(1.3)

cube = result_mL2{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('proton physicalDose with mL objective')
zoom(1.3)

cube = result_mL2{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photon physicalDose with mL objective')
zoom(1.3)

cube = physDose_mL2;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total physicalDose with mL objective')
zoom(1.3)

%% subplot only BED

cube = BED_without_proton;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('proton BED without objective')
zoom(1.3)

cube = BED_without_photon;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photon BED without objective')
zoom(1.3)

cube = BED_without_combi;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total BED without objective')
zoom(1.3)

cube = BED_DD_proton;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('proton BED with DD objective')
zoom(1.3)

cube = BED_DD_photon;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photon BED with DD objective')
zoom(1.3)

cube = BED_DD_combi;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total BED with DD objective')
zoom(1.3)

cube = BED_mL2_proton;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('proton BED with mL objective')
zoom(1.3)

cube = BED_mL2_photon;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photon BED with mL objective')
zoom(1.3)

cube = BED_mL2_combi;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total BED with mL objective')
zoom(1.3)


%% with DD in OAR2
clear
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
cst{4,2}    = 'Margin';
cst{4,3}    = 'OAR';
cst{4,4}{1} = find(mVOIEnlarged);
cst{4,5}    = cst{1,5};

cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40)); 
cst{4,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,2));


% changing alphaX
cst{2,5}.alphaX = 0.5;
%%
% changing alphaX
%cst{2,5}.alphaX = 0.1;

cst{3,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(800,60));
cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% meta information for treatment plan (1) 
pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
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
pln(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
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
pln(2).numOfFractions  = 5;
pln(2).radiationMode   = 'carbon';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [45 0 -45]; % [?] ;
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
pln(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
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

%pln = pln(1);

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
%resultGUI = matRad_calcResultGUIstruct(result);
% cst_secondD = cst;
% dij_secondD = dij2D;
% result_secondD = result2D;
% resultGUI_secondD = result2D{1,1};
result_DD = result;
%
save("3_CarbonProton_Beam_MarginDD_TG119_MixedModalities_withDD.mat","-v7.3")

%%
clear
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
cst{4,2}    = 'Margin';
cst{4,3}    = 'OAR';
cst{4,4}{1} = find(mVOIEnlarged);
cst{4,5}    = cst{1,5};

cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40)); 
cst{4,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,3));

%%
% changing alphaX
cst{2,5}.alphaX = 0.5;

% changing alphaX
%cst{2,5}.alphaX = 0.5;

%cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));
cst{3,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(800,60));
cst{1,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,2));
%cst_third = cst;

% meta information for treatment plan (1) 
pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'carbon';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
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
pln(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
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
pln(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
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

%pln = pln(1);

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
dij = matRad_calcmLETDose(dij,pln);
result = matRad_fluenceOptimizationJO(dij,cst,plnJO);
%resultGUI = matRad_calcResultGUIstruct(result);
% dij_third = dij;
% resultGUI_third = result3{1,1};
result_mL = result;

save("3_CarbonPhoton_Beam_MarginmL_TG119_MixedModalities_withmL.mat","-v7.3")
%%

cube = result_mL{1,1}.dirtyDose;
plane = 3;
slice = 80;
doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('dirty dose with mL objective')
zoom(1.3)

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
 %% DVH Calculation

% physicalDose protons
dvh_without_p = matRad_calcDVH(cst_without,result_without_com{1,1}.physicalDose);
dvh_DD_p = matRad_calcDVH(cst_DD,result_DD{1,1}.physicalDose);
dvh_mL_p = matRad_calcDVH(cst_mL,result_mL{1,1}.physicalDose);

% dirtyDose protons
dvh_without_d = matRad_calcDVH(cst_without,result_without_com{1,1}.dirtyDose);
dvh_DD_d = matRad_calcDVH(cst_DD,result_DD{1,1}.dirtyDose);
dvh_mL_d = matRad_calcDVH(cst_mL,result_mL{1,1}.dirtyDose);

% LET protons
dvh_without_L = matRad_calcDVH(cst_without,result_without_com{1,1}.LET);
dvh_DD_L = matRad_calcDVH(cst_DD,result_DD{1,1}.LET);
dvh_mL_L = matRad_calcDVH(cst_mL,result_mL{1,1}.LET);

% effect
dvh_without_e = matRad_calcDVH(cst_without,result_without_com{1,1}.effect);
dvh_DD_e = matRad_calcDVH(cst_DD,result_DD{1,1}.effect);
dvh_mL_e = matRad_calcDVH(cst_mL,result_mL{1,1}.effect);
%% Total for effect Optimization
% physicalDose
% dvh_without_e = 
% dirtyDose
% LET
%% Total for physicalDose Optimization
% physicalDose
dvh_Pwithout_p = matRad_calcDVH(cst_without,physDose_without);
dvh_PDD_p = matRad_calcDVH(cst_DD,physDose_DD);
dvh_PmL2_p = matRad_calcDVH(cst_mL,physDose_mL2);
% dirtyDose
dvh_Pwithout_d = matRad_calcDVH(cst_without,dirtyDose_without);
dvh_PDD_d = matRad_calcDVH(cst_DD,dirtyDose_DD);
dvh_PmL2_d = matRad_calcDVH(cst_mL,dirtyDose_mL2);
% LET
dvh_Pwithout_L = matRad_calcDVH(cst_without,LET_without);
dvh_PDD_L = matRad_calcDVH(cst_DD,LET_DD);
dvh_PmL2_L = matRad_calcDVH(cst_mL,LET_mL2);

%% indicator wrapper
% physicalDose only proton
qi_Prowithout_p  = matRad_calcQualityIndicators(cst_without,pln(1),result_without{1}.physicalDose,[],[]);
qi_ProDD_p  = matRad_calcQualityIndicators(cst_DD,pln(1),result_DD{1}.physicalDose,[],[]);
qi_PromL2_p  = matRad_calcQualityIndicators(cst_mL,pln(1),result_mL2{1}.physicalDose,[],[]);

% dirtyDose only proton
qi_Prowithout_d  = matRad_calcQualityIndicators(cst_without,pln(1),result_without{1}.dirtyDose,[],[]);
qi_ProDD_d  = matRad_calcQualityIndicators(cst_DD,pln(1),result_DD{1}.dirtyDose,[],[]);
qi_PromL2_d  = matRad_calcQualityIndicators(cst_mL,pln(1),result_mL2{1}.dirtyDose,[],[]);

% LET only proton
qi_Prowithout_L  = matRad_calcQualityIndicators(cst_without,pln(1),result_without{1}.LET,[],[]);
qi_ProDD_L  = matRad_calcQualityIndicators(cst_DD,pln(1),result_DD{1}.LET,[],[]);
qi_PromL2_L  = matRad_calcQualityIndicators(cst_mL,pln(1),result_mL2{1}.LET,[],[]);

% physicalDose only photon
qi_Phowithout_p  = matRad_calcQualityIndicators(cst_without,pln(2),result_without{2}.physicalDose,[],[]);
qi_PhoDD_p  = matRad_calcQualityIndicators(cst_DD,pln(2),result_DD{2}.physicalDose,[],[]);
qi_PhomL2_p  = matRad_calcQualityIndicators(cst_mL,pln(2),result_mL2{2}.physicalDose,[],[]);

% total physicalDose
qi_without_p  = matRad_calcQualityIndicators(cst_without,pln(3),physDose_without,[],[]);
qi_DD_p  = matRad_calcQualityIndicators(cst_DD,pln(3),physDose_DD,[],[]);
qi_mL2_p  = matRad_calcQualityIndicators(cst_mL,pln(3),physDose_mL2,[],[]);

% dirtyDose only proton
qi_without_d  = matRad_calcQualityIndicators(cst_without,pln(3),dirtyDose_without,[],[]);
qi_DD_d  = matRad_calcQualityIndicators(cst_DD,pln(3),dirtyDose_DD,[],[]);
qi_mL2_d  = matRad_calcQualityIndicators(cst_mL,pln(3),dirtyDose_mL2,[],[]);

% LET only proton
qi_without_L  = matRad_calcQualityIndicators(cst_without,pln(3),LET_without,[],[]);
qi_DD_L  = matRad_calcQualityIndicators(cst_DD,pln(3),LET_DD,[],[]);
qi_mL2_L  = matRad_calcQualityIndicators(cst_mL,pln(3),LET_mL2,[],[]);

%% plotting
color = ['r', 'b', 'm' , 'g'];
figure
subplot(3,1,1)
for i = [2 3 4]
    plot(dvh_Pwithout_p(i).doseGrid,dvh_Pwithout_p(i).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',1.5)
    hold on
    plot(dvh_PDD_p(i).doseGrid,dvh_PDD_p(i).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1.5)
    hold on
    plot(dvh_PmL_p(i).doseGrid,dvh_PmL_p(i).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',1.5)
end
title('DVH: blue:Target, magenta:Body, green:Core_Big')
xlabel('physical dose')
ylabel('Volume in %')
legend('without','DD','mL')

subplot(3,1,2)
for i = [2 3 4]
    plot(dvh_Pwithout_d(i).doseGrid,dvh_Pwithout_d(i).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',1.5)
    hold on
    plot(dvh_PDD_d(i).doseGrid,dvh_PDD_d(i).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1.5)
    hold on
    plot(dvh_PmL_d(i).doseGrid,dvh_PmL_d(i).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',1.5)
end
xlabel('dirty dose in Gy')
ylabel('Volume in %')
legend('without','DD','mL')

subplot(3,1,3)
for i = [2 3 4]
    plot(dvh_Pwithout_L(i).doseGrid,dvh_Pwithout_L(i).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',1.5)
    hold on
    plot(dvh_PDD_L(i).doseGrid,dvh_PDD_L(i).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1.5)
    hold on
    plot(dvh_PmL_L(i).doseGrid,dvh_PmL_L(i).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',1.5)
end
xlabel('LET in keV/µm')
ylabel('Volume in %')
legend('without','DD','mL')

% subplot(4,1,4)
% for i = [2 3 4]
%     plot(dvh_without_e(i).doseGrid,dvh_without_e(i).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',1.5)
%     hold on
%     plot(dvh_DD_e(i).doseGrid,dvh_DD_e(i).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1.5)
%     hold on
%     plot(dvh_mL_e(i).doseGrid,dvh_mL_e(i).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',1.5)
% end
% xlabel('effect')
% ylabel('Volume in %')
% legend('without','DD','mL')
%%
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