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


%%

% meta information for treatment plan (1) 
pln_1(1).numOfFractions  = 5;
pln_1(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_1(1).machine         = 'Generic';

% beam geometry settings
pln_1(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_1(1).propStf.gantryAngles    = 0; %[-45 45];%[45 0 -45]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_1(1).propStf.couchAngles     = zeros(numel(pln_1(1).propStf.gantryAngles),1); % [?] ; 
pln_1(1).propStf.numOfBeams      = numel(pln_1(1).propStf.gantryAngles);
pln_1(1).propStf.isoCenter       = ones(pln_1(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln_1(1).propDoseCalc.calcLET = 1;

pln_1(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_1(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_1(1).propOpt.spatioTemp      = 0;
pln_1(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_1(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_1(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_1(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_1(1).bioParam = matRad_bioModel(pln_1(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_1(1).multScen = matRad_multScen(ct,scenGenType);

%pln = pln(1);
%%
% meta information for treatment plan (2) 
pln_1(2).numOfFractions  = 25;
pln_1(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_1(2).machine         = 'Generic';

% beam geometry settings
pln_1(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_1(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_1(2).propStf.couchAngles     = zeros(numel(pln_1(2).propStf.gantryAngles),1);  % [?] ; 
pln_1(2).propStf.numOfBeams      = numel(pln_1(2).propStf.gantryAngles);
pln_1(2).propStf.isoCenter       = ones(pln_1(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln_1(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_1(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_1(2).propOpt.spatioTemp      = 0;
pln_1(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_1(2).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln_1(2).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln_1(2).propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_1(2).bioParam = matRad_bioModel(pln_1(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_1(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);
%%
% Plan Wrapper
%plnJO_1beam_test = matRad_plnWrapper(pln_1);
% Stf Wrapper
%stf_1beamDD100_test = matRad_stfWrapper(ct,cst,plnJO_1beam_test);
stf_LxD100_1beam = matRad_generateStf(ct,cst,pln_1(1));

%%
dij_LxD1 = matRad_calcParticleDose(ct,stf_LxD100_1beam,pln_1,cst);
dij_LxD1 = matRad_calcDirtyDose(2,dij_LxD1,pln_1);
%dij_LxD1 = matRad_calcmLETDose(dij_LxD1,pln_1);
%% Dij Calculation
dij_1beam_L60 = matRad_calcCombiDose(ct,stf_1beam_L60,plnJO_1beam_test,cst,false);
% Dirty Dose Calculation
dij_1beam_L60 = matRad_calcDirtyDose(2,dij_1beam_L60,pln_1);
dij_1beam_L60 = matRad_calcmLETDose(dij_1beam_L60,pln_1);
dij_1beam_L60.precon = 0;
%%
%[result_LxD1_1beam,optimizer] = matRad_fluenceOptimizationJO(dij_1beam_L60,cst,plnJO_1beam_test);
[result_LxD1_1beam,optimizer] = matRad_fluenceOptimization(dij_LxD1,cst,pln_1);

%%
% DD 500 1 beam
physDose_DD500_1 = result_overDD500_1beam{1,1}.physicalDose * 5 + result_overDD500_1beam{1,2}.physicalDose * 25;
RBExD_DD500_1 = result_overDD500_1beam{1,1}.physicalDose * 1.1 * 5 + result_overDD500_1beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_DD500_1 = result_overDD500_1beam{1,1}.dirtyDose * 5 + result_overDD500_1beam{1,2}.dirtyDose * 25;
LET_DD500_1 = (result_overDD500_1beam{1,1}.LET .* result_overDD500_1beam{1,1}.physicalDose * 5 + 0.3 * result_overDD500_1beam{1,2}.physicalDose * 25)./physDose_DD500_1;

% DD 500 2 beams
physDose_DD500_2 = result_overDD500_2beam{1,1}.physicalDose * 5 + result_overDD500_2beam{1,2}.physicalDose * 25;
RBExD_DD500_2 = result_overDD500_2beam{1,1}.physicalDose * 1.1 * 5 + result_overDD500_2beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_DD500_2 = result_overDD500_2beam{1,1}.dirtyDose * 5 + result_overDD500_2beam{1,2}.dirtyDose * 25;
LET_DD500_2 = (result_overDD500_2beam{1,1}.LET .* result_overDD500_2beam{1,1}.physicalDose * 5 + 0.3 * result_overDD500_2beam{1,2}.physicalDose * 25) ./ physDose_DD500_2;

% DD 6 1 beam
physDose_DD6_1 = result_overDD6_1beam{1,1}.physicalDose * 5 + result_overDD6_1beam{1,2}.physicalDose * 25;
RBExD_DD6_1 = result_overDD6_1beam{1,1}.physicalDose * 1.1 * 5 + result_overDD6_1beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_DD6_1 = result_overDD6_1beam{1,1}.dirtyDose * 5 + result_overDD6_1beam{1,2}.dirtyDose * 25;
LET_DD6_1 = (result_overDD6_1beam{1,1}.LET .* result_overDD6_1beam{1,1}.physicalDose * 5 + 0.3 * result_overDD6_1beam{1,2}.physicalDose * 25)./physDose_DD6_1;

% DD 6 2 beams
physDose_DD6_2 = result_overDD6_2beam{1,1}.physicalDose * 5 + result_overDD6_2beam{1,2}.physicalDose * 25;
RBExD_DD6_2 = result_overDD6_2beam{1,1}.physicalDose * 1.1 * 5 + result_overDD6_2beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_DD6_2 = result_overDD6_2beam{1,1}.dirtyDose * 5 + result_overDD6_2beam{1,2}.dirtyDose * 25;
LET_DD6_2 = (result_overDD6_2beam{1,1}.LET .* result_overDD6_2beam{1,1}.physicalDose * 5 + 0.3 * result_overDD6_2beam{1,2}.physicalDose * 25)./physDose_DD6_2;

% DD 100 1 beam
physDose_DD100_1 = result_overDD_1beam{1,1}.physicalDose * 5 + result_overDD_1beam{1,2}.physicalDose * 25;
RBExD_DD100_1 = result_overDD_1beam{1,1}.physicalDose * 1.1 * 5 + result_overDD_1beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_DD100_1 = result_overDD_1beam{1,1}.dirtyDose * 5 + result_overDD_1beam{1,2}.dirtyDose * 25;
LET_DD100_1 = (result_overDD_1beam{1,1}.LET .* result_overDD_1beam{1,1}.physicalDose * 5 + 0.3 * result_overDD_1beam{1,2}.physicalDose * 25)./physDose_DD100_1;

% DD 100 2 beams
physDose_DD100_2 = result_overDD_2beams{1,1}.physicalDose * 5 + result_overDD_2beams{1,2}.physicalDose * 25;
RBExD_DD100_2 = result_overDD_2beams{1,1}.physicalDose * 1.1 * 5 + result_overDD_2beams{1,2}.physicalDose * 1.1 * 25;
dirtyDose_DD100_2 = result_overDD_2beams{1,1}.dirtyDose * 5 + result_overDD_2beams{1,2}.dirtyDose * 25;
LET_DD100_2 = (result_overDD_2beams{1,1}.LET .* result_overDD_2beams{1,1}.physicalDose * 5 + 0.3 * result_overDD_2beams{1,2}.physicalDose * 25)./physDose_DD100_2;

% LxD 3 1 beam 
physDose_LxD3_1 = result_overLxD3_1beam{1,1}.physicalDose * 5 + result_overLxD3_1beam{1,2}.physicalDose * 25;
RBExD_LxD3_1 = result_overLxD3_1beam{1,1}.physicalDose * 1.1 * 5 + result_overLxD3_1beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_LxD3_1 = result_overLxD3_1beam{1,1}.dirtyDose * 5 + result_overLxD3_1beam{1,2}.dirtyDose * 25;
LET_LxD3_1 = (result_overLxD3_1beam{1,1}.LET .* result_overLxD3_1beam{1,1}.physicalDose * 5 + 0.3 * result_overLxD3_1beam{1,2}.physicalDose * 25)./physDose_LxD3_1;

% LxD 3 2 beams
physDose_LxD3_2 = result_overLxD3_2beam{1,1}.physicalDose * 5 + result_overLxD3_2beam{1,2}.physicalDose * 25;
RBExD_LxD3_2 = result_overLxD3_2beam{1,1}.physicalDose * 1.1 * 5 + result_overLxD3_2beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_LxD3_2 = result_overLxD3_2beam{1,1}.dirtyDose * 5 + result_overLxD3_2beam{1,2}.dirtyDose * 25;
LET_LxD3_2 = (result_overLxD3_2beam{1,1}.LET .* result_overLxD3_2beam{1,1}.physicalDose * 5 + 0.3 * result_overLxD3_2beam{1,2}.physicalDose * 25)./physDose_LxD3_2;

% LxD 60 1 beam
physDose_LxD60_1 = result_overLxD60_1beam{1,1}.physicalDose * 5 + result_overLxD60_1beam{1,2}.physicalDose * 25;
RBExD_LxD60_1 = result_overLxD60_1beam{1,1}.physicalDose * 1.1 * 5 + result_overLxD60_1beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_LxD60_1 = result_overLxD60_1beam{1,1}.dirtyDose * 5 + result_overLxD60_1beam{1,2}.dirtyDose * 25;
LET_LxD60_1 = (result_overLxD60_1beam{1,1}.LET .* result_overLxD60_1beam{1,1}.physicalDose * 5 + 0.3 * result_overLxD60_1beam{1,2}.physicalDose * 25)./physDose_LxD60_1;

% LxD 60 2 beams
physDose_LxD60_2 = result_overLxD60_2beams{1,1}.physicalDose * 5 + result_overLxD60_2beams{1,2}.physicalDose * 25;
RBExD_LxD60_2 = result_overLxD60_2beams{1,1}.physicalDose * 1.1 * 5 + result_overLxD60_2beams{1,2}.physicalDose * 1.1 * 25;
dirtyDose_LxD60_2 = result_overLxD60_2beams{1,1}.dirtyDose * 5 + result_overLxD60_2beams{1,2}.dirtyDose * 25;
LET_LxD60_2 = (result_overLxD60_2beams{1,1}.LET .* result_overLxD60_2beams{1,1}.physicalDose * 5 + 0.3 * result_overLxD60_2beams{1,2}.physicalDose * 25)./physDose_LxD60_2;

% LxD 100 1 beam
physDose_LxD100_1 = result_overLxD100_1beam{1,1}.physicalDose * 5 + result_overLxD100_1beam{1,2}.physicalDose * 25;
RBExD_LxD100_1 = result_overLxD100_1beam{1,1}.physicalDose * 1.1 * 5 + result_overLxD100_1beam{1,2}.physicalDose * 1.1 * 25;
dirtyDose_LxD100_1 = result_overLxD100_1beam{1,1}.dirtyDose * 5 + result_overLxD100_1beam{1,2}.dirtyDose * 25;
LET_LxD100_1 = (result_overLxD100_1beam{1,1}.LET .* result_overLxD100_1beam{1,1}.physicalDose * 5 + 0.3 * result_overLxD100_1beam{1,2}.physicalDose * 25)./physDose_LxD100_1;

% LxD 100 2 beams
physDose_LxD100_2 = result_overLxD100_2beams{1,1}.physicalDose * 5 + result_overLxD100_2beams{1,2}.physicalDose * 25;
RBExD_LxD100_2 = result_overLxD100_2beams{1,1}.physicalDose * 1.1 * 5 + result_overLxD100_2beams{1,2}.physicalDose * 1.1 * 25;
dirtyDose_LxD100_2 = result_overLxD100_2beams{1,1}.dirtyDose * 5 + result_overLxD100_2beams{1,2}.dirtyDose * 25;
LET_LxD100_2 = (result_overLxD100_2beams{1,1}.LET .* result_overLxD100_2beams{1,1}.physicalDose * 5 + 0.3 * result_overLxD100_2beams{1,2}.physicalDose * 25)./physDose_LxD100_2;

%% Creating images old 3 beams
cube1 = result_overDD{1,1}.physicalDose * 5;
cube2 = result_overDD{1,2}.physicalDose *25;
cube3 = physDose_O100;
cube4 = result_overDD{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overDD{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_O100;
cube7 = result_overDD{1,1}.dirtyDose * 5;
cube8 = result_overDD{1,2}.dirtyDose * 25;
cube9 = dirtyDose_O100;
% cube10 = (result_overDD{1,1}.LET .* result_overDD{1,1}.physicalDose * 5)/result_overDD{1,1}.physicalDose;
% cube11 = 0.3 * result_overDD{1,2}.physicalDose * 25 /result_overDD{1,2}.physicalDose;
% cube12 = LET_O100;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)

%% Creating images 0°
% DD 6 1 beam
cube1 = result_overDD6_1beam{1,1}.physicalDose * 5;
cube2 = result_overDD6_1beam{1,2}.physicalDose *25;
cube3 = physDose_DD6_1;
cube4 = result_overDD6_1beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overDD6_1beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_DD6_1;
cube7 = result_overDD6_1beam{1,1}.dirtyDose * 5;
cube8 = result_overDD6_1beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_DD6_1;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons DD6 1 beam 0°')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)

% DD 6 2 beams
cube1 = result_overDD6_2beam{1,1}.physicalDose * 5;
cube2 = result_overDD6_2beam{1,2}.physicalDose *25;
cube3 = physDose_DD6_2;
cube4 = result_overDD6_2beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overDD6_2beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_DD6_2;
cube7 = result_overDD6_2beam{1,1}.dirtyDose * 5;
cube8 = result_overDD6_2beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_DD6_2;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons DD6 2 beams [-45° 45°]')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)


% DD 100 1 beam
cube1 = result_overDD_1beam{1,1}.physicalDose * 5;
cube2 = result_overDD_1beam{1,2}.physicalDose *25;
cube3 = physDose_DD100_1;
cube4 = result_overDD_1beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overDD_1beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_DD100_1;
cube7 = result_overDD_1beam{1,1}.dirtyDose * 5;
cube8 = result_overDD_1beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_DD100_1;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons DD100 1 beam 0°')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)

% DD 100 2 beams
cube1 = result_overDD_2beams{1,1}.physicalDose * 5;
cube2 = result_overDD_2beams{1,2}.physicalDose *25;
cube3 = physDose_DD100_2;
cube4 = result_overDD_2beams{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overDD_2beams{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_DD100_2;
cube7 = result_overDD_2beams{1,1}.dirtyDose * 5;
cube8 = result_overDD_2beams{1,2}.dirtyDose * 25;
cube9 = dirtyDose_DD100_2;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons DD100 2 beams [-45° 45°]')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)


% DD 500 1 beam
cube1 = result_overDD500_1beam{1,1}.physicalDose * 5;
cube2 = result_overDD500_1beam{1,2}.physicalDose *25;
cube3 = physDose_DD500_1;
cube4 = result_overDD500_1beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overDD500_1beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_DD500_1;
cube7 = result_overDD500_1beam{1,1}.dirtyDose * 5;
cube8 = result_overDD500_1beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_DD500_1;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons DD500 1 beam 0°')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)

% DD 500 2 beams
cube1 = result_overDD500_2beam{1,1}.physicalDose * 5;
cube2 = result_overDD500_2beam{1,2}.physicalDose *25;
cube3 = physDose_DD500_2;
cube4 = result_overDD500_2beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overDD500_2beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_DD500_2;
cube7 = result_overDD500_2beam{1,1}.dirtyDose * 5;
cube8 = result_overDD500_2beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_DD500_2;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons DD500 2 beam [-45° 45°]')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)


% LxD 3 1 beam
cube1 = result_overLxD3_1beam{1,1}.physicalDose * 5;
cube2 = result_overLxD3_1beam{1,2}.physicalDose *25;
cube3 = physDose_LxD3_1;
cube4 = result_overLxD3_1beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overLxD3_1beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_LxD3_1;
cube7 = result_overLxD3_1beam{1,1}.dirtyDose * 5;
cube8 = result_overLxD3_1beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_LxD3_1;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons LxD3 1 beam 0°')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)

% LxD 3 2 beams
cube1 = result_overLxD3_2beam{1,1}.physicalDose * 5;
cube2 = result_overLxD3_2beam{1,2}.physicalDose *25;
cube3 = physDose_LxD3_2;
cube4 = result_overLxD3_2beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overLxD3_2beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_LxD3_2;
cube7 = result_overLxD3_2beam{1,1}.dirtyDose * 5;
cube8 = result_overLxD3_2beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_LxD3_2;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons LxD3 2 beams [-45° 45°]')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)


% LxD 60 1 beam
cube1 = result_overLxD60_1beam{1,1}.physicalDose * 5;
cube2 = result_overLxD60_1beam{1,2}.physicalDose *25;
cube3 = physDose_LxD60_1;
cube4 = result_overLxD60_1beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overLxD60_1beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_LxD60_1;
cube7 = result_overLxD60_1beam{1,1}.dirtyDose * 5;
cube8 = result_overLxD60_1beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_LxD60_1;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons LxD60 1 beam 0°')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)

% LxD 60 2 beams
cube1 = result_overLxD60_2beams{1,1}.physicalDose * 5;
cube2 = result_overLxD60_2beams{1,2}.physicalDose *25;
cube3 = physDose_LxD60_2;
cube4 = result_overLxD60_2beams{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overLxD60_2beams{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_LxD60_2;
cube7 = result_overLxD60_2beams{1,1}.dirtyDose * 5;
cube8 = result_overLxD60_2beams{1,2}.dirtyDose * 25;
cube9 = dirtyDose_LxD60_2;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons LxD60 2 beams [-45° 45°]')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)


% LxD 100 1 beam
cube1 = result_overLxD100_1beam{1,1}.physicalDose * 5;
cube2 = result_overLxD100_1beam{1,2}.physicalDose *25;
cube3 = physDose_LxD100_1;
cube4 = result_overLxD100_1beam{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overLxD100_1beam{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_LxD100_1;
cube7 = result_overLxD100_1beam{1,1}.dirtyDose * 5;
cube8 = result_overLxD100_1beam{1,2}.dirtyDose * 25;
cube9 = dirtyDose_LxD100_1;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons LxD100 1 beam 0°')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)

% LxD 100 2 beams
cube1 = result_overLxD100_2beams{1,1}.physicalDose * 5;
cube2 = result_overLxD100_2beams{1,2}.physicalDose *25;
cube3 = physDose_LxD100_2;
cube4 = result_overLxD100_2beams{1,1}.physicalDose * 1.1 * 5;
cube5 = result_overLxD100_2beams{1,2}.physicalDose * 1.1 *25;
cube6 = RBExD_LxD100_2;
cube7 = result_overLxD100_2beams{1,1}.dirtyDose * 5;
cube8 = result_overLxD100_2beams{1,2}.dirtyDose * 25;
cube9 = dirtyDose_LxD100_2;

plane = 3;
slice = 80;
doseWindow = [0 70];

figure,
subplot(3,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('protons LxD100 2 beams [-45° 45°]')
zoom(2)
subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photons')
zoom(2)
subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total')
zoom(2)
subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)
subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube9,plane,slice,[],[],colorcube,[],doseWindow,[]);
zoom(2)


%%
cube = result_without{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [0 70];

subplot(3,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Reference proton dose')
zoom(4)

cube = result_without{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [0 70];

subplot(3,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Reference photon dose')
zoom(4)

cube = physDose_U100;
plane = 3;
slice = 80;
doseWindow = [0 70];

subplot(3,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total dose DD U(100,30)')
zoom(4)

cube = result_U100{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [0 70];

subplot(3,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('proton dose DD U(100,30)')
zoom(4)

cube = result_U100{1,1}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [0 70];

subplot(3,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photon dose DD U(100,30)')
zoom(4)

cube = physDose_U100O100;
plane = 3;
slice = 80;
doseWindow = [0 70];

subplot(3,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total dose DD U(100,30), O(100,0)')
zoom(4)

cube = result_U100O100{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [0 70];

subplot(3,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('proton dose U(100,30), O(100,0)')
zoom(4)

cube = result_U100O100{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [0 70];

subplot(3,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('photon dose U(100,30), O(100,0)')
zoom(4)
%%
cube = physDose_DD10;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(5,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis DD(10,0)')
zoom(4)

% diffProton30 = result_without{1,1}.physicalDose - result_DD30{1,1}.physicalDose;
% diffProton10 = result_without{1,1}.physicalDose - result_DD10{1,1}.physicalDose;
% diffPhoton30 = result_without{1,2}.physicalDose - result_DD30{1,2}.physicalDose;
% diffPhoton10 = result_without{1,2}.physicalDose - result_DD10{1,2}.physicalDose;
% diff30Proton10 = result_DD30{1,1}.physicalDose - result_DD10{1,1}.physicalDose;
% diff30Photon10 = result_DD30{1,2}.physicalDose - result_DD10{1,2}.physicalDose;


cube = result_DD30{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [0 60];

subplot(5,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis DD(30,0)')
zoom(4)

cube = result_DD30{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [0 60];

subplot(5,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis DD(30,0)')
zoom(4)

cube = physDose_DD30;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(5,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis DD(30,0)')
zoom(4)

cube = result_DD{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [0 60];

subplot(5,3,10)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis DD(100,0)')
zoom(4)

cube = result_DD{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [0 60];

subplot(5,3,11)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis DD(100,0)')
zoom(4)

cube = physDose_DD;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(5,3,12)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis DD(100,0)')
zoom(4)

cube = result_DD300{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [0 60];

subplot(5,3,13)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis DD(300,0)')
zoom(4)

cube = result_DD300{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [0 60];

subplot(5,3,14)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis DD(300,0)')
zoom(4)

cube = physDose_DD300;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(5,3,15)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis DD(300,0)')
zoom(4)
%%
array = [physDose_mL6,physDose_mL10,physDose_mL30,physDose_mL];
cube = result_without{1,1}.physicalDose *5;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

figure,
subplot(5,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenzprotonendosis')
zoom(4)

cube = result_without{1,2}.physicalDose *25;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenzphotonendosis')
zoom(4)

cube = physDose_without;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(array(:))];

subplot(5,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz-Gesamtdosis')
zoom(4)

cube = result_mL6{1,1}.physicalDose *5;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis LxD(6,0)')
zoom(4)

cube = result_mL6{1,2}.physicalDose *25;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis LxD(6,0)')
zoom(4)

cube = physDose_mL6;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(array(:))];

subplot(5,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis LxD(6,0)')
zoom(4)

cube = result_mL10{1,1}.physicalDose *5;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis LxD(10,0)')
zoom(4)

cube = result_mL10{1,2}.physicalDose *25;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis LxD(10,0)')
zoom(4)

cube = physDose_mL10;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(array(:))];

subplot(5,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis LxD(10,0)')
zoom(4)

cube = result_mL30{1,1}.physicalDose *5;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,10)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis LxD(30,0)')
zoom(4)

cube = result_mL30{1,2}.physicalDose *25;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,11)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis LxD(30,0)')
zoom(4)

cube = physDose_mL30;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(array(:))];

subplot(5,3,12)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis LxD(30,0)')
zoom(4)

cube = result_mL{1,1}.physicalDose *5;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,13)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis LxD(100,0)')
zoom(4)

cube = result_mL{1,2}.physicalDose *25;
plane = 3;
slice = 80;
doseWindow = [0 max(array(:))];

subplot(5,3,14)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis LxD(100,0)')
zoom(4)

cube = physDose_mL;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(array(:))];

subplot(5,3,15)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis LxD(100,0)')
zoom(4)

%%
cube = result_without{1,1}.physicalDose ;
plane = 3;
slice = 80;
doseWindow = [0 12];

figure,
subplot(7,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenzprotonendosis')
zoom(4)

cube = result_without{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 3];

subplot(7,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenzphotonendosis')
zoom(4)

cube = physDose_without;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(7,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz-Gesamtdosis')
zoom(4)

cube = result_DD10{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 12];

subplot(7,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis DD(10,0)')
zoom(4)

cube = result_DD10{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 3];

subplot(7,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis DD(10,0)')
zoom(4)

cube = physDose_DD10;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(7,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis DD(10,0)')
zoom(4)

cube = result_mL10{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 12];

subplot(7,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis LxD(10,0)')
zoom(4)

cube = result_mL10{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 3];

subplot(7,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis LxD(10,0)')
zoom(4)

cube = physDose_mL10;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(7,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis LxD(10,0)')
zoom(4)

cube = result_DD30{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 12];

subplot(7,3,10)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis DD(30,0)')
zoom(4)

cube = result_DD30{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 3];

subplot(7,3,11)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis DD(30,0)')
zoom(4)

cube = physDose_DD30;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(7,3,12)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis DD(30,0)')
zoom(4)

cube = result_mL30{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 12];

subplot(7,3,13)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis LxD(30,0)')
zoom(4)

cube = result_mL30{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 3];

subplot(7,3,14)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis LxD(30,0)')
zoom(4)

cube = physDose_mL30;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(7,3,15)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis LxD(30,0)')
zoom(4)

cube = result_DD{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 12];

subplot(7,3,16)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis DD(100,0)')
zoom(4)

cube = result_DD{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 3];

subplot(7,3,17)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis DD(100,0)')
zoom(4)

cube = physDose_DD;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(7,3,18)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis DD(100,0)')
zoom(4)

cube = result_mL{1,1}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 12];

subplot(7,3,19)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis LxD(100,0)')
zoom(4)

cube = result_mL{1,2}.physicalDose;
plane = 3;
slice = 80;
doseWindow = [0 3];

subplot(7,3,20)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis LxD(100,0)')
zoom(4)

cube = physDose_mL;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];

subplot(7,3,21)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis LxD(100,0)')
zoom(4)



%% difference map total
diff = physDose_DD10 - physDose_DD30;
cube = diff;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];
figure
matRad_plotSliceWrapper(gca,ct,cst_DD10,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('total physDose DD10 - DD30')

%% BED
% BED = D*(1+d/(alpha/beta))
BED_without_proton = result_without{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_without{1,1}.physicalDose./(0.1./0.05));
BED_DD_proton = result_DD{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_DD{1,1}.physicalDose./(0.1./0.05));
BED_DD10_proton = result_DD10{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_DD10{1,1}.physicalDose./(0.1./0.05));
BED_DD30_proton = result_DD30{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_DD30{1,1}.physicalDose./(0.1./0.05));
BED_DD300_proton = result_DD300{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_DD300{1,1}.physicalDose./(0.1./0.05));
BED_mL_proton = result_mL{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_mL{1,1}.physicalDose./(0.1./0.05));
BED_mL6_proton = result_mL6{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_mL6{1,1}.physicalDose./(0.1./0.05));
BED_mL10_proton = result_mL10{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_mL10{1,1}.physicalDose./(0.1./0.05));
BED_mL30_proton = result_mL30{1,1}.physicalDose * pln_1(1).numOfFractions .*(1 + result_mL30{1,1}.physicalDose./(0.1./0.05));

BED_without_photon = result_without{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_without{1,2}.physicalDose./(0.1./0.05));
BED_DD_photon = result_DD{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_DD{1,2}.physicalDose./(0.1./0.05));
BED_DD10_photon = result_DD10{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_DD10{1,2}.physicalDose./(0.1./0.05));
BED_DD30_photon = result_DD30{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_DD30{1,2}.physicalDose./(0.1./0.05));
BED_DD300_photon = result_DD300{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_DD300{1,2}.physicalDose./(0.1./0.05));
BED_mL_photon = result_mL{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_mL{1,2}.physicalDose./(0.1./0.05));
BED_mL6_photon = result_mL6{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_mL6{1,2}.physicalDose./(0.1./0.05));
BED_mL10_photon = result_mL10{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_mL10{1,2}.physicalDose./(0.1./0.05));
BED_mL30_photon = result_mL30{1,2}.physicalDose * pln_1(2).numOfFractions .*(1 + result_mL30{1,2}.physicalDose./(0.1./0.05));

BED_without_combi = BED_without_proton + BED_without_photon;
BED_DD_combi = BED_DD_proton + BED_DD_photon;
BED_DD10_combi = BED_DD10_proton + BED_DD10_photon;
BED_DD30_combi = BED_DD30_proton + BED_DD30_photon;
BED_DD300_combi = BED_DD300_proton + BED_DD300_photon;
BED_mL_combi = BED_mL_proton + BED_mL_photon;
BED_mL6_combi = BED_mL6_proton + BED_mL6_photon;
BED_mL10_combi = BED_mL10_proton + BED_mL10_photon;
BED_mL30_combi = BED_mL30_proton + BED_mL30_photon;

%% subplot total
array = [physDose_DD,physDose_DD10,physDose_DD30,physDose_mL,physDose_mL10,physDose_mL30,physDose_without];
array_D = [dirtyDose_DD,dirtyDose_DD30,dirtyDose_DD10,dirtyDose_mL,dirtyDose_mL10,dirtyDose_mL30,dirtyDose_without];
array_L = [LET_DD,LET_DD30,LET_DD10,LET_mL,LET_mL10,LET_mL30,LET_without];

cube = physDose_without;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

figure,
subplot(7,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz-Gesamtdosis')
zoom(4)

cube = dirtyDose_without;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(7,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz-Gesamt-DirtyDose')
zoom(4)

cube = LET_without;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(7,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz-Gesamt-LET')
zoom(4)

%array_D = [dirtyDose_without,dirtyDose_DD,dirtyDose_mL];

cube = physDose_DD10;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(7,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit DD(10,0)')
zoom(4)

cube = dirtyDose_DD10;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(7,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit DD(10,0)')
zoom(4)

cube = LET_DD10;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(7,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit DD(10,0)')
zoom(4)

cube = physDose_mL10;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(7,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(10,0)')
zoom(4)

cube = dirtyDose_mL10;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(7,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit LxD(10,0)')
zoom(4)

cube = LET_mL10;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(7,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit LxD(10,0)')
zoom(4)

cube = physDose_DD30;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(7,3,10)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit DD(30,0)')
zoom(4)

cube = dirtyDose_DD30;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(7,3,11)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit DD(30,0)')
zoom(4)

cube = LET_DD30;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(7,3,12)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit DD(30,0)')
zoom(4)

cube = physDose_mL30;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(7,3,13)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(30,0)')
zoom(4)

cube = dirtyDose_mL30;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(7,3,14)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit LxD(30,0)')
zoom(4)

cube = LET_mL30;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(7,3,15)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit LxD(30,0)')
zoom(4)

cube = physDose_DD;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(7,3,16)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit DD(100,0)')
zoom(4)

cube = dirtyDose_DD;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(7,3,17)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit DD(100,0)')
zoom(4)

cube = LET_DD;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(7,3,18)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit DD(100,0)')
zoom(4)

cube = physDose_mL;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(7,3,19)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(100,0)')
zoom(4)

cube = dirtyDose_mL;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(7,3,20)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit LxD(100,0)')
zoom(4)

cube = LET_mL;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(7,3,21)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit LxD(100,0)')
zoom(4)

%%
%array_L = [LET_without_2,LET_DD_2,LET_mL_2];

cube = physDose_DD300;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(5,3,13)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit DD(300,0)')
zoom(4)

cube = dirtyDose_DD300;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(5,3,14)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit DD(300,0)')
zoom(4)

cube = LET_DD300;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(5,3,15)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit DD(300,0)')
zoom(4)
%%
array = [physDose_mL,physDose_mL10,physDose_mL30,physDose_mL6,physDose_without];
array_D = [dirtyDose_mL,dirtyDose_mL30,dirtyDose_mL10,dirtyDose_mL6,dirtyDose_without];
array_L = [LET_mL,LET_mL10,LET_mL30,LET_mL6,LET_without];

cube = physDose_without;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

figure,
subplot(5,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz-Gesamtdosis')
zoom(4)

cube = dirtyDose_without;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(5,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz-Gesamt-DirtyDose')
zoom(4)

cube = LET_without;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(5,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz-Gesamt-LET')
zoom(4)

cube = physDose_mL6;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(5,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(6,0)')
zoom(4)

cube = dirtyDose_mL6;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(5,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit LxD(6,0)')
zoom(4)

cube = LET_mL6;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(5,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit LxD(6,0)')
zoom(4)

cube = physDose_mL10;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(5,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(10,0)')
zoom(4)

cube = dirtyDose_mL10;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(5,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit LxD(10,0)')
zoom(4)

cube = LET_mL10;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(5,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit LxD(10,0)')
zoom(4)

cube = physDose_mL30;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(5,3,10)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(30,0)')
zoom(4)

cube = dirtyDose_mL30;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(5,3,11)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit LxD(30,0)')
zoom(4)

cube = LET_mL30;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(5,3,12)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit LxD(30,0)')
zoom(4)

cube = physDose_mL;
plane = 3;
slice = 80;
doseWindow = [min(array(:)) max(array(:))];

subplot(5,3,13)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(100,0)')
zoom(4)

cube = dirtyDose_mL;
plane = 3;
slice = 80;
doseWindow = [min(array_D(:)) max(array_D(:))];

subplot(5,3,14)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-DirtyDose mit LxD(100,0)')
zoom(4)

cube = LET_mL;
plane = 3;
slice = 80;
doseWindow = [min(array_L(:)) max(array_L(:))];

subplot(5,3,15)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-LET mit LxD(100,0)')
zoom(4)

%% subplot only physicalDose
array = [result_without{1,1}.physicalDose*5,result_without{1,2}.physicalDose*25];
array_p = [result_DD{1,1}.physicalDose*5,result_DD{1,2}.physicalDose*25];
array_p2 = [result_DD300{1,1}.physicalDose*5,result_DD300{1,2}.physicalDose*25];
array_p3 = [result_DD30{1,1}.physicalDose*5,result_DD30{1,2}.physicalDose*25];
array_L = [result_mL{1,1}.physicalDose*5,result_mL{1,2}.physicalDose*25];
array_L2 = [result_mL6{1,1}.physicalDose*5,result_mL6{1,2}.physicalDose*25];
array_t = [physDose_without,physDose_DD,physDose_DD300,physDose_mL,physDose_mL6];
array_g = [array,array_p,array_p2,array_p3,array_L,array_L2,array_t];
array_g2 = [result_DD{1,1}.physicalDose*5,result_DD300{1,1}.physicalDose*5,result_DD30{1,1}.physicalDose*5,result_mL{1,1}.physicalDose*5,result_mL6{1,1}.physicalDose*5];

cube = result_without{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

figure,
subplot(6,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz Protonendosis')
zoom(4)

cube = result_without{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz Photonendosis')
zoom(4)

cube = physDose_without;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz Gesamtdosis')
zoom(4)

cube = result_DD{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [min(array_g2(:)) max(array_g2(:))];

subplot(6,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis mit DD(100,0)')
zoom(4)

cube = result_DD{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis mit DD(100,0)')
zoom(4)

cube = physDose_DD;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit DD(100,0)')
zoom(4)

cube = result_DD300{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [min(array_g2(:)) max(array_g2(:))];

subplot(6,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis mit DD(300,0)')
zoom(4)

cube = result_DD300{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis mit DD(300,0)')
zoom(4)

cube = physDose_DD300;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit DD(300,0)')
zoom(4)

cube = result_DD30{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [min(array_g2(:)) max(array_g2(:))];

subplot(6,3,10)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis mit DD(30,0)')
zoom(4)

cube = result_DD30{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,11)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis mit DD(30,0)')
zoom(4)

cube = physDose_DD30;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,12)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit DD(30,0)')
zoom(4)

cube = result_mL{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [min(array_g2(:)) max(array_g2(:))];

subplot(6,3,13)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis mit LxD(100,0)')
zoom(4)

cube = result_mL{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,14)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis mit LxD(100,0)')
zoom(4)

cube = physDose_mL;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,15)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(100,0)')
zoom(4)

cube = result_mL6{1,1}.physicalDose * 5;
plane = 3;
slice = 80;
doseWindow = [min(array_g2(:)) max(array_g2(:))];

subplot(6,3,16)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonendosis mit LxD(6,0)')
zoom(4)

cube = result_mL6{1,2}.physicalDose * 25;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,17)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonendosis mit LxD(6,0)')
zoom(4)

cube = physDose_mL6;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(6,3,18)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamtdosis mit LxD(6,0)')
zoom(4)


%% subplot only BED
array = [BED_without_proton,BED_without_photon];
array_d = [BED_mL_proton,BED_mL_photon];
array_d2 = [BED_mL10_proton,BED_mL10_photon];
array_d3 = [BED_mL30_proton,BED_mL30_photon];
array_d4 = [BED_DD_proton,BED_DD_photon];
array_d5 = [BED_DD10_proton,BED_DD10_photon];
array_d6 = [BED_DD30_proton,BED_DD30_photon];
array_d7 = [BED_without_combi,BED_mL_combi,BED_mL10_combi,BED_mL30_combi,BED_DD_combi,BED_DD10_combi,BED_DD30_combi];
array_g = [array,array_d,array_d2,array_d3,array_d4,array_d5,array_d6,array_d7];
% array_L = [BED_mL_proton,BED_mL_photon];
% array_L2 = [BED_mL6_proton,BED_mL6_photon];
% array_t = [BED_without_combi,BED_DD_combi,BED_mL_combi];
% array_g = [array,array_d,array_d2,array_L,array_L2,array_t];

cube = BED_without_proton;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

figure,
subplot(7,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz Proton BED')
zoom(4)

cube = BED_without_photon;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz Photon BED')
zoom(4)

cube = BED_without_combi;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Referenz Gesamt-BED')
zoom(4)

cube = BED_DD10_proton;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonen BED mit DD(10,0)')
zoom(4)

cube = BED_DD10_photon;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonen BED mit DD(10,0)')
zoom(4)

cube = BED_DD10_combi;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-BED mit DD(10,0)')
zoom(4)

cube = BED_mL10_proton;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonen BED mit LxD(10,0)')
zoom(4)

cube = BED_mL10_photon;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonen BED mit LxD(10,0)')
zoom(4)

cube = BED_mL10_combi;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,9)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-BED mit LxD(10,0)')
zoom(4)

cube = BED_DD30_proton;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,10)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonen BED mit DD(30,0)')
zoom(4)

cube = BED_DD30_photon;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,11)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonen BED mit DD(30,0)')
zoom(4)

cube = BED_DD30_combi;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,12)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-BED mit DD(30,0)')
zoom(4)

cube = BED_mL30_proton;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,13)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonen BED mit LxD(30,0)')
zoom(4)

cube = BED_mL30_photon;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,14)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonen BED mit LxD(30,0)')
zoom(4)

cube = BED_mL30_combi;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,15)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-BED mit LxD(30,0)')
zoom(4)

cube = BED_DD_proton;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,16)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonen BED mit DD(100,0)')
zoom(4)

cube = BED_DD_photon;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,17)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonen BED mit DD(100,0)')
zoom(4)

cube = BED_DD_combi;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,18)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-BED mit DD(100,0)')
zoom(4)

cube = BED_mL_proton;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,19)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Protonen BED mit LxD(100,0)')
zoom(4)

cube = BED_mL_photon;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,20)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Photonen BED mit LxD(100,0)')
zoom(4)

cube = BED_mL_combi;
plane = 3;
slice = 80;
doseWindow = [min(array_g(:)) max(array_g(:))];

subplot(7,3,21)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Gesamt-BED mit LxD(100,0)')
zoom(4)

%% Total for physicalDose Optimization
% physicalDose
dvh_Pwithout_p = matRad_calcDVH(cst_without,physDose_without);
dvh_PDD_p = matRad_calcDVH(cst_DD,physDose_DD);
%dvh_P3DD_p = matRad_calcDVH(cst_DD,physDose_3DD);
dvh_PmL_p = matRad_calcDVH(cst_mL,physDose_mL);
dvh_PDD300_p = matRad_calcDVH(cst_DD300,physDose_DD300);
dvh_PDD30_p = matRad_calcDVH(cst_DD30,physDose_DD30);
dvh_PDD10_p = matRad_calcDVH(cst_DD10,physDose_DD10);
dvh_PmL6_p = matRad_calcDVH(cst_mL6,physDose_mL6);
dvh_PmL10_p = matRad_calcDVH(cst_mL10,physDose_mL10);
dvh_PmL30_p = matRad_calcDVH(cst_mL30,physDose_mL30);
%dvh_PmL2_p = matRad_calcDVH(cst_mL2,physDose_mL2);
% dirtyDose
dvh_Pwithout_d = matRad_calcDVH(cst_without,dirtyDose_without);
dvh_PDD_d = matRad_calcDVH(cst_DD,dirtyDose_DD);
% dvh_P3DD_d = matRad_calcDVH(cst_DD,dirtyDose_3DD);
dvh_PmL_d = matRad_calcDVH(cst_mL,dirtyDose_mL);
dvh_PDD300_d = matRad_calcDVH(cst_DD300,dirtyDose_DD300);
dvh_PDD30_d = matRad_calcDVH(cst_DD30,dirtyDose_DD30);
dvh_PDD10_d = matRad_calcDVH(cst_DD10,dirtyDose_DD10);
dvh_PmL6_d = matRad_calcDVH(cst_mL6,dirtyDose_mL6);
dvh_PmL10_d = matRad_calcDVH(cst_mL10,dirtyDose_mL10);
dvh_PmL30_d = matRad_calcDVH(cst_mL30,dirtyDose_mL30);
% dvh_PmL2_d = matRad_calcDVH(cst_mL2,dirtyDose_mL2);
% LET
dvh_Pwithout_L = matRad_calcDVH(cst_without,LET_without);
dvh_PDD_L = matRad_calcDVH(cst_DD,LET_DD);
% dvh_P3DD_L = matRad_calcDVH(cst_DD,LET_3DD);
dvh_PmL_L = matRad_calcDVH(cst_mL,LET_mL);
dvh_PDD300_L = matRad_calcDVH(cst_DD300,LET_DD300);
dvh_PDD30_L = matRad_calcDVH(cst_DD30,LET_DD30);
dvh_PDD10_L = matRad_calcDVH(cst_DD10,LET_DD10);
dvh_PmL6_L = matRad_calcDVH(cst_mL6,LET_mL6);
dvh_PmL10_L = matRad_calcDVH(cst_mL10,LET_mL10);
dvh_PmL30_L = matRad_calcDVH(cst_mL30,LET_mL30);
% dvh_PmL2_L = matRad_calcDVH(cst_mL2,LET_mL2);



%% indicator wrapper
% physicalDose only proton
qi_Prowithout_p  = matRad_calcQualityIndicators(cst_without,pln_1(1),result_without{1}.physicalDose,[],[]);
qi_ProDD_p  = matRad_calcQualityIndicators(cst_DD,pln_1(1),result_DD{1}.physicalDose,[],[]);
qi_PromL_p  = matRad_calcQualityIndicators(cst_mL,pln_1(1),result_mL{1}.physicalDose,[],[]);
qi_ProDD300_p  = matRad_calcQualityIndicators(cst_DD300,pln_1(1),result_DD300{1}.physicalDose,[],[]);
qi_ProDD30_p  = matRad_calcQualityIndicators(cst_DD30,pln_1(1),result_DD30{1}.physicalDose,[],[]);
qi_ProDD10_p = matRad_calcQualityIndicators(cst_DD10,pln_1(1),result_DD10{1}.physicalDose,[],[]);
qi_PromL6_p  = matRad_calcQualityIndicators(cst_mL6,pln_1(1),result_mL6{1}.physicalDose,[],[]);
qi_PromL10_p  = matRad_calcQualityIndicators(cst_mL10,pln_1(1),result_mL10{1}.physicalDose,[],[]);
qi_PromL30_p  = matRad_calcQualityIndicators(cst_mL30,pln_1(1),result_mL30{1}.physicalDose,[],[]);

% dirtyDose only proton
qi_Prowithout_d  = matRad_calcQualityIndicators(cst_without,pln_o(1),result_without{1}.dirtyDose,[],[]);
qi_ProDD_d  = matRad_calcQualityIndicators(cst_DD,pln_o(1),result_DD{1}.dirtyDose,[],[]);
qi_PromL_d  = matRad_calcQualityIndicators(cst_mL,pln_o(1),result_mL{1}.dirtyDose,[],[]);
qi_ProDD300_d  = matRad_calcQualityIndicators(cst_DD300,pln_o(1),result_DD300{1}.dirtyDose,[],[]);
qi_ProDD30_d  = matRad_calcQualityIndicators(cst_DD30,pln_o(1),result_DD30{1}.dirtyDose,[],[]);
qi_ProDD10_d  = matRad_calcQualityIndicators(cst_DD10,pln_o(1),result_DD10{1}.dirtyDose,[],[]);
qi_PromL6_d  = matRad_calcQualityIndicators(cst_mL6,pln_o(1),result_mL6{1}.dirtyDose,[],[]);
qi_PromL10_p  = matRad_calcQualityIndicators(cst_mL10,pln_1(1),result_mL10{1}.physicalDose,[],[]);
qi_PromL30_p  = matRad_calcQualityIndicators(cst_mL30,pln_1(1),result_mL30{1}.physicalDose,[],[]);

% LET only proton
qi_Prowithout_L  = matRad_calcQualityIndicators(cst_without,pln_o(1),result_without{1}.LET,[],[]);
qi_ProDD_L  = matRad_calcQualityIndicators(cst_DD,pln_o(1),result_DD{1}.LET,[],[]);
qi_PromL_L  = matRad_calcQualityIndicators(cst_mL,pln_o(1),result_mL{1}.LET,[],[]);
qi_ProDD300_L  = matRad_calcQualityIndicators(cst_DD300,pln_o(1),result_DD300{1}.LET,[],[]);
qi_ProDD30_L  = matRad_calcQualityIndicators(cst_DD30,pln_o(1),result_DD30{1}.LET,[],[]);
qi_ProDD10_L  = matRad_calcQualityIndicators(cst_DD10,pln_o(1),result_DD10{1}.LET,[],[]);
qi_PromL6_L  = matRad_calcQualityIndicators(cst_mL6,pln_o(1),result_mL6{1}.LET,[],[]);
qi_PromL10_p  = matRad_calcQualityIndicators(cst_mL10,pln_1(1),result_mL10{1}.physicalDose,[],[]);
qi_PromL30_p  = matRad_calcQualityIndicators(cst_mL30,pln_1(1),result_mL30{1}.physicalDose,[],[]);


% physicalDose only photon
qi_Phowithout_p  = matRad_calcQualityIndicators(cst_without,pln_o(2),result_without{2}.physicalDose,[],[]);
qi_PhoDD_p  = matRad_calcQualityIndicators(cst_DD,pln_o(2),result_DD{2}.physicalDose,[],[]);
qi_PhomL_p  = matRad_calcQualityIndicators(cst_mL,pln_o(2),result_mL{2}.physicalDose,[],[]);
qi_PhoDD300_p  = matRad_calcQualityIndicators(cst_DD300,pln_o(2),result_DD300{2}.physicalDose,[],[]);
qi_PhoDD30_p  = matRad_calcQualityIndicators(cst_DD30,pln_o(2),result_DD30{2}.physicalDose,[],[]);
qi_PhoDD10_p  = matRad_calcQualityIndicators(cst_DD10,pln_o(2),result_DD10{2}.physicalDose,[],[]);
qi_PhomL6_p  = matRad_calcQualityIndicators(cst_mL6,pln_o(2),result_mL6{2}.physicalDose,[],[]);
%%
% total physicalDose
qi_without_p  = matRad_calcQualityIndicators(cst_without,pln3(3),physDose_without,[],[]);
qi_DD_p  = matRad_calcQualityIndicators(cst_DD,pln3(3),physDose_DD,[],[]);
qi_DD300_p  = matRad_calcQualityIndicators(cst_DD300,pln3(3),physDose_DD300,[],[]);
qi_DD30_p  = matRad_calcQualityIndicators(cst_DD30,pln3(3),physDose_DD30,[],[]);
qi_DD10_p  = matRad_calcQualityIndicators(cst_DD10,pln3(3),physDose_DD10,[],[]);
qi_mL_p  = matRad_calcQualityIndicators(cst_mL,pln3(3),physDose_mL,[],[]);
qi_mL6_p  = matRad_calcQualityIndicators(cst_mL6,pln3(3),physDose_mL6,[],[]);
qi_mL10_p  = matRad_calcQualityIndicators(cst_mL10,pln3(3),physDose_mL10,[],[]);
qi_mL30_p  = matRad_calcQualityIndicators(cst_mL30,pln3(3),physDose_mL30,[],[]);


% dirtyDose only proton
qi_without_d  = matRad_calcQualityIndicators(cst_without,pln3(3),dirtyDose_without,[],[]);
qi_DD_d  = matRad_calcQualityIndicators(cst_DD,pln3(3),dirtyDose_DD,[],[]);
qi_DD300_d  = matRad_calcQualityIndicators(cst_DD300,pln3(3),dirtyDose_DD300,[],[]);
qi_DD30_d  = matRad_calcQualityIndicators(cst_DD30,pln3(3),dirtyDose_DD30,[],[]);
qi_DD10_d  = matRad_calcQualityIndicators(cst_DD10,pln3(3),dirtyDose_DD10,[],[]);
qi_mL_d  = matRad_calcQualityIndicators(cst_mL,pln3(3),dirtyDose_mL,[],[]);
qi_mL6_d  = matRad_calcQualityIndicators(cst_mL6,pln3(3),dirtyDose_mL6,[],[]);
qi_mL10_d  = matRad_calcQualityIndicators(cst_mL10,pln3(3),dirtyDose_mL10,[],[]);
qi_mL30_d  = matRad_calcQualityIndicators(cst_mL30,pln3(3),dirtyDose_mL30,[],[]);

% LET only proton
qi_without_L  = matRad_calcQualityIndicators(cst_without,pln3(3),LET_without,[],[]);
qi_DD_L  = matRad_calcQualityIndicators(cst_DD,pln3(3),LET_DD,[],[]);
qi_DD300_L  = matRad_calcQualityIndicators(cst_DD300,pln3(3),LET_DD300,[],[]);
qi_DD30_L  = matRad_calcQualityIndicators(cst_DD30,pln3(3),LET_DD30,[],[]);
qi_DD10_L  = matRad_calcQualityIndicators(cst_DD10,pln3(3),LET_DD10,[],[]);
qi_mL_L  = matRad_calcQualityIndicators(cst_mL,pln3(3),LET_mL,[],[]);
qi_mL6_L  = matRad_calcQualityIndicators(cst_mL6,pln3(3),LET_mL6,[],[]);
qi_mL10_L  = matRad_calcQualityIndicators(cst_mL10,pln3(3),LET_mL10,[],[]);
qi_mL30_L  = matRad_calcQualityIndicators(cst_mL30,pln3(1),LET_mL30,[],[]);

%% plotting
color = ['b', 'm', 'g' , 'r'];
figure
subplot(3,1,1)
structure = [2 3 4];
for i = [1 2 3]
    plot(dvh_Pwithout_p(structure(i)).doseGrid,dvh_Pwithout_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',2)
    hold on
    plot(dvh_PDD10_p(structure(i)).doseGrid,dvh_PDD10_p(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',2.5)
    hold on
    plot(dvh_PDD30_p(structure(i)).doseGrid,dvh_PDD30_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',2)
    hold on
    plot(dvh_PDD_p(structure(i)).doseGrid,dvh_PDD_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',2)
    % hold on
    % plot(dvh_PDD300_p(structure(i)).doseGrid,dvh_PDD300_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1)
    hold on
    %plot(dvh_PmL6_p(structure(i)).doseGrid,dvh_PmL6_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',1)
    % hold on
    plot(dvh_PmL10_p(structure(i)).doseGrid,dvh_PmL10_p(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',1.5)
    hold on
    plot(dvh_PmL30_p(structure(i)).doseGrid,dvh_PmL30_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1)
    hold on
    plot(dvh_PmL_p(structure(i)).doseGrid,dvh_PmL_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',1.5)
end
title('DVH: blau:PTV, magenta:Body, grün:Core-Big')
xlabel('physikalische Dosis in Gy')
ylabel('Volumen in %')
legend('Referenzplan','DD(10,0)','DD(30,0)','DD(100,0)','LxD(10,0)','LxD(30,0)','LxD(100,0)')%'DD(10,0)','DD(30,0)','DD(100,0)','DD(300,0)')%,'LxD(6,0)','LxD(10,0)','LxD(30,0)','LxD(100,0)')
%legend('Referenzplan','DD(30,0)','DD(10,0)')
subplot(3,1,2)
for i = [1 2 3]
    plot(dvh_Pwithout_d(structure(i)).doseGrid,dvh_Pwithout_d(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',2)
    hold on
    plot(dvh_PDD10_d(structure(i)).doseGrid,dvh_PDD10_d(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',2.5)
    hold on
    plot(dvh_PDD30_d(structure(i)).doseGrid,dvh_PDD30_d(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',2)
    hold on
    plot(dvh_PDD_d(structure(i)).doseGrid,dvh_PDD_d(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',2)
    hold on
    % plot(dvh_PDD300_d(structure(i)).doseGrid,dvh_PDD300_d(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1)
    % hold on
    % plot(dvh_PmL6_d(structure(i)).doseGrid,dvh_PmL6_d(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',1)
    % hold on
    plot(dvh_PmL10_d(structure(i)).doseGrid,dvh_PmL10_d(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',1.5)
    hold on
    plot(dvh_PmL30_d(structure(i)).doseGrid,dvh_PmL30_d(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1)
    hold on
    plot(dvh_PmL_d(structure(i)).doseGrid,dvh_PmL_d(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',1.5)
end
xlabel('DirtyDose in Gy')
ylabel('Volumen in %')
legend('Referenzplan','DD(10,0)','DD(30,0)','DD(100,0)','LxD(10,0)','LxD(30,0)','LxD(100,0)')%,'DD(10,0)','DD(30,0)','DD(100,0)','DD(300,0)')%,'LxD(6,0)','LxD(10,0)','LxD(30,0)','LxD(100,0)')
%legend('Referenzplan','DD(30,0)','DD(10,0)')
subplot(3,1,3)
for i = [1 2 3]
    plot(dvh_Pwithout_L(structure(i)).doseGrid,dvh_Pwithout_L(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',2)
    hold on
    plot(dvh_PDD10_L(structure(i)).doseGrid,dvh_PDD10_L(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',2.5)
    hold on
    plot(dvh_PDD30_L(structure(i)).doseGrid,dvh_PDD30_L(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',2)
    hold on
    plot(dvh_PDD_L(structure(i)).doseGrid,dvh_PDD_L(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',2)
    hold on
    % plot(dvh_PDD300_L(structure(i)).doseGrid,dvh_PDD300_L(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1)
    % hold on
    % plot(dvh_PmL6_L(structure(i)).doseGrid,dvh_PmL6_L(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',1)
    % hold on
    plot(dvh_PmL10_L(structure(i)).doseGrid,dvh_PmL10_L(structure(i)).volumePoints,'Color',color(i),'LineStyle',':','LineWidth',1.5)
    hold on
    plot(dvh_PmL30_L(structure(i)).doseGrid,dvh_PmL30_L(structure(i)).volumePoints,'Color',color(i),'LineStyle','--','LineWidth',1)
    hold on
    plot(dvh_PmL_L(structure(i)).doseGrid,dvh_PmL_L(structure(i)).volumePoints,'Color',color(i),'LineStyle','-.','LineWidth',1.5)
end
xlabel('LET in keV/µm')
ylabel('Volumen in %')
legend('Referenzplan','DD(10,0)','DD(30,0)','DD(100,0)','LxD(10,0)','LxD(30,0)','LxD(100,0)')%'DD(10,0)','DD(30,0)','DD(100,0)','DD(300,0)')%,'LxD(6,0)','LxD(10,0)','LxD(30,0)','LxD(100,0)')
%legend('Referenzplan','DD(30,0)','DD(10,0)')
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
figure
xaxis = [10];
xaxis2 = [30];
xaxis3 = 100;
xaxis4 = 300;
xaxis5 = 6;
yaxis = [qi_DD10_p(2).D_95];
yaxis_b = [qi_DD10_p(4).mean];
yaxis_c = [qi_DD10_d(4).mean];
yaxis_d = [qi_DD10_p(3).mean];

yaxis2 = [qi_DD30_p(2).D_95];
yaxis2_b = [qi_DD30_p(4).mean];
yaxis2_c = [qi_DD30_d(4).mean];
yaxis2_d = [qi_DD30_p(3).mean];

yaxis3 = [qi_DD_p(2).D_95];
yaxis3_b = [qi_DD_p(4).mean];
yaxis3_c = [qi_DD_d(4).mean];
yaxis3_d = [qi_DD_p(3).mean];

yaxis4 = [qi_DD300_p(2).D_95];
yaxis4_b = [qi_DD300_p(4).mean];
yaxis4_c = [qi_DD300_d(4).mean];
yaxis4_d = [qi_DD300_p(3).mean];

yaxis5 = [qi_mL6_p(2).D_95];
yaxis5_b = [qi_mL6_p(4).mean];
yaxis5_c = [qi_mL6_d(4).mean];
yaxis5_d = [qi_mL6_p(3).mean];

yaxis3b = [qi_mL_p(2).D_95];
yaxis3b_b = [qi_mL_p(4).mean];
yaxis3b_c = [qi_mL_d(4).mean];
yaxis3b_d = [qi_mL_p(3).mean];

scatter(xaxis,yaxis,'magenta')
hold on
scatter(xaxis,yaxis_b,'blue')
hold on
scatter(xaxis,yaxis_c,'green')
hold on
scatter(xaxis,yaxis_d,'yellow')
hold on
scatter(xaxis2,yaxis2,'magenta')
hold on
scatter(xaxis2,yaxis2_b,'blue')
hold on
scatter(xaxis2,yaxis2_c,'green')
hold on
scatter(xaxis2,yaxis2_d,'yellow')
hold on
scatter(xaxis3,yaxis3,'magenta')
hold on
scatter(xaxis3,yaxis3_b,'blue')
hold on
scatter(xaxis3,yaxis3_c,'green')
hold on
scatter(xaxis3,yaxis3_d,'yellow')
hold on
scatter(xaxis4,yaxis4,'magenta')
hold on
scatter(xaxis4,yaxis4_b,'blue')
hold on
scatter(xaxis4,yaxis4_c,'green')
hold on
scatter(xaxis4,yaxis4_d,'yellow')

legend('PTV95','OAR dose mean','OAR DD mean','Body dose mean')