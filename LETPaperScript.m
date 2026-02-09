%% My processing script 
matRad_rc;
% clear
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 500000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-06;
load TG119.mat
% load('D:\postDoc\LET_Paper\FinalResult.mat')
% load('D:\postDoc\LET_Paper\OnlyProton.mat')
 %%
%% add core
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

%% add Target margin
cube = zeros(ct.cubeDim);
cube(cst{1,4}{1}) = 1;
vResolution = ct.resolution;
vMargin = [];
vMargin.x = 5;
vMargin.y = 5;
vMargin.z = 5;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin,1);

cst{4,1}    = 3;
cst{4,2}    = 'Core-big';
cst{4,3}    = 'OAR';
cst{4,4}{1} = find(mVOIEnlarged);
cst{4,5}    = cst{1,5};
cst{4,6}    = cst{1,6};

%cst{2,5}.alphaX  = 0.5;
%cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,60));
%% TG119
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));
cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));

cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0)); 
cst{4,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0)); 

cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,60));

cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));
cst{3,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));

%%
cst{3,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,20)); 
cst{13,6}{2} = struct(mLETDoseObjectives.matRad_SquaredUnderdosingmLETDose(100,90)); 
cst{4,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0)); 
cst{13,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,18)); 
% cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,40)); 
 cst{3,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));

 %% Glioma case different alpha/beta ratios

% a/ß for brainstem = 2.1 Gy, we can keep it

% a/ß for glioma (GTV and CTV) = 6 Gy, a = 0.25 Gy^-1, ß = 0.04167 Gy^-2
cst{13,5}.alphaX = 0.25;
cst{14,5}.alphaX = 0.25;
cst{13,5}.betaX = 0.04167;
cst{14,5}.betaX = 0.04167;

%% recreating ref joint plan

cst{3,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));

%%

% % meta information for treatment plan (1) 
pln_onlyPro(1).numOfFractions  = 5;
pln_onlyPro(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln_onlyPro(1).machine         = 'Generic';

% beam geometry settings
pln_onlyPro(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_onlyPro(1).propStf.gantryAngles    = [265 285]; %[90 180 270];% [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln_onlyPro(1).propStf.couchAngles     = zeros(numel(pln_onlyPro(1).propStf.gantryAngles),1); % [?] ; 
pln_onlyPro(1).propStf.numOfBeams      = numel(pln_onlyPro(1).propStf.gantryAngles);
pln_onlyPro(1).propStf.isoCenter       = ones(pln_onlyPro(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln_onlyPro(1).propDoseCalc.calcLET = 1;

pln_onlyPro(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_onlyPro(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_onlyPro(1).propOpt.spatioTemp      = 0;
pln_onlyPro(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_onlyPro(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_onlyPro(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_onlyPro(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
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
pln_onlyPro(1).bioParam = matRad_bioModel(pln_onlyPro(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_onlyPro(1).multScen = matRad_multScen(ct,scenGenType);

%pln = pln(1);

% meta information for treatment plan (2) 
pln_onlyPro(2).numOfFractions  = 25;
pln_onlyPro(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln_onlyPro(2).machine         = 'Generic';

% beam geometry settings
pln_onlyPro(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln_onlyPro(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln_onlyPro(2).propStf.couchAngles     = zeros(numel(pln_onlyPro(2).propStf.gantryAngles),1);  % [?] ; 
pln_onlyPro(2).propStf.numOfBeams      = numel(pln_onlyPro(2).propStf.gantryAngles);
pln_onlyPro(2).propStf.isoCenter       = ones(pln_onlyPro(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln_onlyPro(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln_onlyPro(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln_onlyPro(2).propOpt.spatioTemp      = 0;
pln_onlyPro(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln_onlyPro(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln_onlyPro(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln_onlyPro(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

%quantityOpt  = ['effect'];     % options: physicalDose, effect, RBExD
quantityOpt = 'physicalDose';
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln_onlyPro(2).bioParam = matRad_bioModel(pln_onlyPro(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln_onlyPro(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);

%pln = pln(2);

% Plan Wrapper
plnJO = matRad_plnWrapper(pln_onlyPro);
%% Stf Wrapper
stf = matRad_stfWrapper(ct,cst,plnJO);

%% Dij Calculation
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
% Dirty Dose Calculation
dij = matRad_calcDirtyDose(2,dij,pln_onlyPro);
dij = matRad_calcmLETDose(dij,pln_onlyPro);

%% compress optimized results into one cube 
%
dij.precon = 1;
% dij.wInit = [result_pre{1}.w; result_pre{2}.w];
[result_GlioDDu,optimizer] = matRad_fluenceOptimizationJO(dij,cst,plnJO);
%% Visualization
slice = 65;

photon_plan = result_GlioRefJ{2};
proton_plan = result_GlioRefJ{1};

effect = @(n,d,a,b)(n.*(a.*d+b.*d.*d))
rbeEffect = effect(30,2,0.1,0.05)

quantityOpt = 'physicalDose';

totalPlan = pln_onlyPro(1).numOfFractions.*proton_plan.(quantityOpt) + pln_onlyPro(2).numOfFractions.*photon_plan.(quantityOpt);

f = figure;
subplot(1,3,1);
    imagesc(proton_plan.(quantityOpt)(:,:,slice).*pln_onlyPro(1).numOfFractions);
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Proton Plan');
subplot(1,3,2);
    imagesc(photon_plan.(quantityOpt)(:,:,slice).*pln_onlyPro(2).numOfFractions);
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Photon Plan');
subplot(1,3,3);
    imagesc(totalPlan(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Total Plan');


%%
plane = 3;
slice = 65;
cube = result_GlioRefJ{1}.physicalDose .* pln_onlyPro(1).numOfFractions;
doseWindow = [0 70];%[0 max(cube(:))];
%isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
figure,
subplot(1,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow);
title(['proton physical dose'])
zoom(1.5)

subplot(1,3,2)
cube = result_pre{2}.physicalDose .* pln_onlyPro(2).numOfFractions;
doseWindow = [0 70];
%isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow);
title(['photon physical dose'])
zoom(1.5)

totalPlan = pln_onlyPro(1).numOfFractions.*proton_plan.physicalDose + pln_onlyPro(2).numOfFractions.*photon_plan.physicalDose;

subplot(1,3,3)
cube = totalPlan;
doseWindow = [0 70];%[0 max(cube(:))];
%isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow);
title('total physical dose')
zoom(1.5)

%%
if iscell(result_pre)
[resultGUI_DDOU, result_DDOU] = matRad_accumulateCubesMixMod(result_pre, pln_onlyPro,ct);
end
% pln = pln(1) ;
% resultGUI.totalEffect = resultGUI.mod1effect + resultGUI.mod2effect;

%% DVH

%% plotting
color = ['b', 'm', 'g' , 'r'];
figure
subplot(3,1,1)
structure = [2 3 4];
for i = [1 2 3]
    plot(dvh_Pwithout_p(structure(i)).doseGrid,dvh_Pwithout_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',2)
    hold on
end

%% Visualization
slice = 65;
ResultCell = resultLET_Target;
photon_plan = ResultCell{2};
proton_plan = ResultCell{1};
quantityOpt = 'effect';
totalPlanOnly = pln_DDovernunder(1).numOfFractions.*RefJoint{1,1}.(quantityOpt) + pln_DDovernunder(2).numOfFractions.*RefJoint{1,2}.(quantityOpt);
% matRad_calcQualityIndicators(cst,pln,totalPlan)
totalOnlyPlan = pln_DDovernunder(1).numOfFractions.*onlyProton{1,1}.(quantityOpt);


%% %% ficures slice 

plane = 3;
slice = 65;
cube = RefJoint{1,2}.RBExD;
doseWindow = [0 2.3];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
figure,
subplot(4,4,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title(['Ref Plan photon RBExD'])
zoom(1.5)

subplot(4,4,3)
cube = result_DDOU{1,2}.RBExD;
doseWindow = [0 2.3];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title(['joint DD both photon fractions RBExD'])
zoom(1.5)

subplot(4,4,4)
cube = result_LETOU{1,2}.RBExD;
doseWindow = [0 2.3];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('joint LET both photon fractions RBExD')
zoom(1.5)

plane = 3;
slice = 65;
cube = onlyProton{1,1}.RBExD;
doseWindow = [0 6.8];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];

subplot(4,4,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title(['Only Proton Plan RBExD' ])
zoom(1.5)

subplot(4,4,6)
cube = RefJoint{1,1}.RBExD;
doseWindow = [0 6.8];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title(['Ref Joint Proton RBExD'])
zoom(1.5)

subplot(4,4,7)
cube = result_DDOU{1,1}.RBExD;
doseWindow = [0 6.8];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('joint DD both proton fractions RBExD')
zoom(1.5)

subplot(4,4,8)
cube = result_LETOU{1,1}.RBExD;
doseWindow = [0 6.8];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('joint LET both proton fractions RBExD')
zoom(1.5)

subplot(4,4,9)
cube = onlyProtonTotal;
doseWindow = [0 155];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('Proton Plan BED')
zoom(1.5)

subplot(4,4,10)
cube = refJointTotal;
doseWindow = [0 155];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('Ref Plan BED')
zoom(1.5)

subplot(4,4,11)
cube = DDTotalOU;
doseWindow = [0 155];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('DD both BED')
zoom(1.5)

subplot(4,4,12)
cube = LETTotalOU;
doseWindow = [0 155];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('LET both BED')
zoom(1.5)

subplot(4,4,13)
cube = onlyProton{1,1}.dirtyDose;
doseWindow = [0 4];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('Proton only Plan dirtyDose')
zoom(1.5)

subplot(4,4,14)
cube = RefJoint{1,1}.dirtyDose;
doseWindow = [0 4];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('Ref Plan dirtyDose')
zoom(1.5)

subplot(4,4,15)
cube = result_DDOU{1,1}.dirtyDose;
doseWindow = [0 4];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('DD both Plan dirtyDose')
zoom(1.5)

subplot(4,4,16)
cube = result_LETOU{1,1}.dirtyDose;
doseWindow = [0 4];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('LET both Plan dirtyDose')
zoom(1.5)

%%
figure
cube = proton_plan_DDovernunder.dirtyDose;
doseWindow = [0 4.8 ];%[0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('Proton dirty dose both (threshold = 2 keV \mum^{-1})')
zoom(1.5)



%% calc DVH for the current cube 
photonPlan = (30.*result_photon{1,1}.effect)./0.1; % + pln(2).numOfFractions.*result_DDovernunder{1,2}.effect)./0.1;
LETTotalOU = (pln_onlyPro(1).numOfFractions.*result_LETOU{1,1}.effect + pln_onlyPro(2).numOfFractions.*result_LETOU{1,2}.effect)./0.1;
DDTotalOU = (pln_onlyPro(1).numOfFractions.*result_DDOU{1,1}.effect + pln_onlyPro(2).numOfFractions.*result_DDOU{1,2}.effect)./0.1;

LETTotalover = (pln_onlyPro(1).numOfFractions.*result_LET_over{1,1}.effect + pln_onlyPro(2).numOfFractions.*result_LET_over{1,2}.effect)./0.1;
DDTotalover = (pln_onlyPro(1).numOfFractions.*result_DD_over{1,1}.effect + pln_onlyPro(2).numOfFractions.*result_DD_over{1,2}.effect)./0.1;

LETTotalunder = (pln_onlyPro(1).numOfFractions.*resultLET_Target{1,1}.effect + pln_onlyPro(2).numOfFractions.*resultLET_Target{1,2}.effect)./0.1;
DDTotalunder = (pln_onlyPro(1).numOfFractions.*resultDD_Target{1,1}.effect + pln_onlyPro(2).numOfFractions.*resultDD_Target{1,2}.effect)./0.1;

% DDOproton = (pln(1).numOfFractions.*result_protonDDO{1,1}.effect)./0.1;
% LETOproton = (pln(1).numOfFractions.*result_protonLETO{1,1}.effect)./0.1;

onlyProtonTotal = effect(30,onlyProton{1,1}.RBExD,0.1,0.05)./0.1  ;
refJointTotal = (pln_onlyPro(1).numOfFractions.*RefJoint{1,1}.effect + pln_onlyPro(2).numOfFractions.*RefJoint{1,2}.effect)./0.1; 
% resultDVH  = matRad_calcDVH(cst,totalPlan,'cum');
% LET_DVH = resultDVH;
%%
cst{3,5}.visibleColor = [ 0 0.7 0];
cst{13,5}.visibleColor = [ 0.2 0.45 1];

matRad_setOverlapPriorities(cst,ct)
%% total BED DVH
%DDtotalDVH  = matRad_calcDVH(cst,DDTotal,'cum');
%LETtotalDVH = matRad_calcDVH(cst,LETTotal,'cum');
onlyProtonDVH = matRad_calcDVH(cst,onlyProtonTotal,'cum');
refJointDVH = matRad_calcDVH(cst,refJointTotal,'cum');

DDOUDVH  = matRad_calcDVH(cst,DDTotalOU,'cum');
LETOUDVH  = matRad_calcDVH(cst,LETTotalOU,'cum');

DDoverDVH  = matRad_calcDVH(cst,DDTotalover,'cum');
LEToverDVH  = matRad_calcDVH(cst,LETTotalover,'cum');

DDunderDVH  = matRad_calcDVH(cst,DDTotalunder,'cum');
LETunderDVH  = matRad_calcDVH(cst,LETTotalunder,'cum');

%photonPlanDVH = matRad_calcDVH(cst,photonPlan,'cum');
%%
figure,
%subplot(1,2,1)
vois = [3,12,13,14,15];
for i  = vois
    % plot(DDOUDVH(i).doseGrid,DDOUDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    % hold on
    plot(LETOUDVH(i).doseGrid,LETOUDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(LEToverDVH(i).doseGrid,LEToverDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
    % plot(DDunderDVH(i).doseGrid,DDunderDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    % hold on
    plot(LETunderDVH(i).doseGrid,LETunderDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    hold on
    % plot(LETtotalDVH(i).doseGrid,LETtotalDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    % hold on
    % plot(onlyProtonDVH(i).doseGrid,onlyProtonDVH(i).volumePoints,'LineWidth',2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    % hold on
    %plot(photonPlanDVH(i).doseGrid,photonPlanDVH(i).volumePoints,'LineWidth',2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    %hold on
    plot(refJointDVH(i).doseGrid,refJointDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'joint extrema';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'joint SO';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.2);
names{c} = 'joint SU';
c = c+1;
% leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
% names{c} = 'joint LET SU';
% c = c+1;
% leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',2);
% names{c} = 'only Proton';
% c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
names{c} = 'Ref. Joint';
c = c+1;
axis([0,160,0,100])
xlabel('BED [Gy]')
ylabel('Volume [%]')
legend ([leg{:}],names ) 
set(gca,'FontSize',14)
grid on

%% color and line style vice versa
figure;
vois = [3,12,13,14,15];

% Farben nach Plan
planColors = containers.Map();
planColors('DDOUDVH')       = [0.8 0 0];      
planColors('LETOUDVH')      = [0 0.6 0.8];    
planColors('DDoverDVH')     = [0 0.4 0];      
planColors('LEToverDVH')    = [0.5 0 0.5];    
planColors('onlyProtonDVH') = [1 0.5 0];      
planColors('refJointDVH')   = [0 0 0];        

% Darstellung je VOI
voistyles = containers.Map({3,12,13,14}, {'-','--',':','-.'}); % Linienstile
voiMarkers = containers.Map({15}, {'o'});                      % VOI 15 mit Marker

% Plot
for i = vois
    if isKey(voistyles, i)
        ls = voistyles(i);
        mk = 'none';  % kein Marker
    elseif isKey(voiMarkers, i)
        ls = 'none';  % kein Linienstil
        mk = voiMarkers(i);
    end
    
    % MarkerIndizes reduziert Markeranzahl (nicht bei linestyle)
    markerSetting = {};
    if ~strcmp(mk,'none')
        markerSetting = {'Marker', mk, 'MarkerIndices', 1:20:1000};
    end

    % Plots je Plan
    plot(DDOUDVH(i).doseGrid,      DDOUDVH(i).volumePoints,      'Color', planColors('DDOUDVH'),       'LineStyle', ls, markerSetting{:}, 'LineWidth', 1.5); hold on
    plot(LETOUDVH(i).doseGrid,     LETOUDVH(i).volumePoints,     'Color', planColors('LETOUDVH'),      'LineStyle', ls, markerSetting{:}, 'LineWidth', 1.5); hold on
    plot(DDoverDVH(i).doseGrid,    DDoverDVH(i).volumePoints,    'Color', planColors('DDoverDVH'),     'LineStyle', ls, markerSetting{:}, 'LineWidth', 1.5); hold on
    plot(LEToverDVH(i).doseGrid,   LEToverDVH(i).volumePoints,   'Color', planColors('LEToverDVH'),    'LineStyle', ls, markerSetting{:}, 'LineWidth', 1.5); hold on
    plot(onlyProtonDVH(i).doseGrid,onlyProtonDVH(i).volumePoints,'Color', planColors('onlyProtonDVH'), 'LineStyle', ls, markerSetting{:}, 'LineWidth', 1.5); hold on
    plot(refJointDVH(i).doseGrid,  refJointDVH(i).volumePoints,  'Color', planColors('refJointDVH'),   'LineStyle', ls, markerSetting{:}, 'LineWidth', 1.5); hold on
end

% Legende Pläne (Farben)
leg = {};
names = {};
plans = {'DDOUDVH','LETOUDVH','DDoverDVH','LEToverDVH','onlyProtonDVH','refJointDVH'};
planNames = {'joint DD SO SU','joint LET SO SU','joint DD SO','joint LET SO','only Proton','Ref. Joint'};

for p = 1:length(plans)
    leg{p} = plot(NaN, 'Color', planColors(plans{p}), 'LineStyle','-', 'LineWidth',1.5);
    names{p} = planNames{p};
end

% Trenner
leg{end+1} = plot(NaN, NaN, 'Color', 'none');
names{end+1} = '--- Organ structures ---';

% Legende VOIs (Stil oder Marker)
for i = vois
    if isKey(voistyles, i)
        leg{end+1} = plot(NaN, 'Color', [0 0 0], 'LineStyle', voistyles(i), 'LineWidth', 1.5);
    elseif isKey(voiMarkers, i)
        leg{end+1} = plot(NaN, 'Color', [0 0 0], 'LineStyle', 'none', 'Marker', voiMarkers(i), 'LineWidth', 1.5);
    end
    names{end+1} = cst{i,2}; % Organname
end

xlabel('BED [Gy]')
ylabel('Volume [%]')
legend([leg{:}], names, 'Location','best');
set(gca,'FontSize',14)
grid on

%% DVH Dirty Dose 

%% Dirty Dose DVH
% DDtotalDVH  = matRad_calcDVH(cst,result_DD_over{1,1}.dirtyDose,'cum');
% LETtotalDVH = matRad_calcDVH(cst,result_LET_over{1,1}.dirtyDose,'cum');
% onlyProtonDVH = matRad_calcDVH(cst,onlyProton{1}.dirtyDose,'cum');
% refJointDVH = matRad_calcDVH(cst,RefJoint{1,1}.dirtyDose,'cum');
% 
DDOUtotalDDDVH  = matRad_calcDVH(cst,DDOU_DBED,'cum');
LETOUtotalDDDVH = matRad_calcDVH(cst,LETOU_DBED,'cum');
DDovertotalDDDVH  = matRad_calcDVH(cst,DDover_DBED,'cum');
LETovertotalDDDVH = matRad_calcDVH(cst,LETover_DBED,'cum');
DDundertotalDDDVH  = matRad_calcDVH(cst,DDunder_DBED,'cum');
LETundertotalDDDVH = matRad_calcDVH(cst,LETunder_DBED,'cum');
onlyProtonDDDVH = matRad_calcDVH(cst,onlyProton_DBED,'cum');
refJointDDDVH = matRad_calcDVH(cst,RefJoint_DBED,'cum');
%%
figure;
%subplot(1,2,2)

vois = [3,13];
for i  = vois
    plot(DDOUtotalDDDVH(i).doseGrid,DDOUtotalDDDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    % plot(LETOUtotalDDDVH(i).doseGrid,LETOUtotalDDDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    % hold on
    plot(DDovertotalDDDVH(i).doseGrid,DDovertotalDDDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
    plot(DDundertotalDDDVH(i).doseGrid,DDundertotalDDDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    hold on
    % plot(LETundertotalDDDVH(i).doseGrid,LETundertotalDDDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    % hold on
    % plot(onlyProtonDDDVH(i).doseGrid,onlyProtonDDDVH(i).volumePoints,'LineWidth',2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    % hold on
    plot(refJointDDDVH(i).doseGrid,refJointDDDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'joint extrema';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'joint SO';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.2);
names{c} = 'joint SU';
c = c+1;
% leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
% names{c} = 'joint LET SU';
% c = c+1;
% leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',2);
% names{c} = 'only Proton';
% c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
names{c} = 'Ref. Joint';
c = c+1;
axis([0,60,0,100])
legend ([leg{:}],names ) 
xlabel(' Dirty BED [Gy]')
ylabel('Volume [%]')
set(gca,'FontSize',14)
 grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.Box = 'on'; 

%% Chat GPT example

figure;
    subplot(1,2,1);
    vois = [3,12,13,14,15];
    ddvois = [3,13];
    for i  = vois
        % plot(DDOUDVH(i).doseGrid,DDOUDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
        % hold on
        plot(DDunderDVH(i).doseGrid,DDunderDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
        hold on
        % plot(DDunderDVH(i).doseGrid,DDunderDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
        % hold on
        plot(LETunderDVH(i).doseGrid,LETunderDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
        hold on
        plot(onlyProtonDVH(i).doseGrid,onlyProtonDVH(i).volumePoints,'LineWidth',2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
        hold on
        plot(refJointDVH(i).doseGrid,refJointDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
        grid on;
        ax = gca;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.Box = 'on';  % keep frame if you want full box look
    end

    hold on
    c =1;
    leg = {};
    names = {};
    %custom legend
    for i = vois
        leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.6);
        hold on
        names{c} = cst{i,2};
        c = c+1;

    end
    leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.6);
    names{c} = 'joint DD';
    c = c+1;
    leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.6);
    names{c} = 'joint LxD';
    c = c+1;
    leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.6);
    names{c} = 'only Proton';
    c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
    % names{c} = 'joint LET SU';
    % c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',2);
    % names{c} = 'only Proton';
    % c = c+1;
    leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.6);
    names{c} = 'Ref. Joint';
    c = c+1;
    axis([0,160,0,100])
    legend ([leg{:}],names )
    ylabel('Volume [%]')
    xlabel('BED [Gy]')
    set(gca,'FontSize',18)
    
    subplot(2,2,2);
    for i  = ddvois
    % plot(DDOUtotalDDDVH(i).doseGrid,DDOUtotalDDDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    % hold on
    plot(DDundertotalDDDVH(i).doseGrid,DDundertotalDDDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    hold on
    plot(LETundertotalDDDVH(i).doseGrid,LETundertotalDDDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
    hold on
    plot(onlyProtonDDDVH(i).doseGrid,onlyProtonDDDVH(i).volumePoints,'LineWidth',2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(refJointDDDVH(i).doseGrid,refJointDDDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.Box = 'on';  % keep frame if you want full box look
    end
    hold on
    c =1;
    leg = {};
    names = {};
    %custom legend
    for i = ddvois
        leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.6);
        hold on
        names{c} = cst{i,2};
        c = c+1;

    end
    leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.6);
    names{c} = 'joint DD';
    c = c+1;
    leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.6);
    names{c} = 'joint LxD';
    c = c+1;
    leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.6);
    names{c} = 'only Proton';
    c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
    % names{c} = 'joint LET SU';
    % c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',2);
    % names{c} = 'only Proton';
    % c = c+1;
    leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.6);
    names{c} = 'Ref. Joint';
    c = c+1;
    axis([0,60,0,100])
    legend ([leg{:}],names )
    % ax.XColor = 'none';
    ax.YColor = 'none';
    xlabel('Dirty BED [Gy]')
    set(gca,'FontSize',18)

    % subplot(2,2,3);
    % for i  = vois
    % plot(LETOUDVH(i).doseGrid,LETOUDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    % hold on
    % plot(LEToverDVH(i).doseGrid,LEToverDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
    % hold on
    % plot(LETunderDVH(i).doseGrid,LETunderDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    % hold on
    % plot(refJointDVH(i).doseGrid,refJointDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    % grid on;
    % ax = gca;
    % ax.XGrid = 'on';
    % ax.YGrid = 'on';
    % ax.Box = 'on';  % keep frame if you want full box look
    % end
    % hold on
    % c =1;
    % leg = {};
    % names = {};
    % %custom legend
    % for i = vois
    %     leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.6);
    %     hold on
    %     names{c} = cst{i,2};
    %     c = c+1;
    % 
    % end
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.6);
    % names{c} = 'joint MinMax';
    % c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.6);
    % names{c} = 'joint SO';
    % c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.6);
    % names{c} = 'joint SU';
    % c = c+1;
    % % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
    % % names{c} = 'joint LET SU';
    % % c = c+1;
    % % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',2);
    % % names{c} = 'only Proton';
    % % c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.6);
    % names{c} = 'Ref. Joint';
    % c = c+1;
    % axis([0,160,0,100])
    % legend ([leg{:}],names )
    % ylabel(['LxD' ...
    %     'Volume [%]'])
    % xlabel('BED [GY]')
    % set(gca,'FontSize',18)
    % 
    % subplot(2,2,4);
    % for i  = ddvois
    % plot(LETOUtotalDDDVH(i).doseGrid,LETOUtotalDDDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    % hold on
    % plot(LETovertotalDDDVH(i).doseGrid,LETovertotalDDDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
    % hold on
    % plot(LETundertotalDDDVH(i).doseGrid,LETundertotalDDDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    % hold on
    % plot(refJointDDDVH(i).doseGrid,refJointDDDVH(i).volumePoints,'LineWidth',1.6,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    % grid on;
    % ax = gca;
    % ax.XGrid = 'on';
    % ax.YGrid = 'on';
    % ax.Box = 'on';  % keep frame if you want full box look
    % end
    % hold on
    % c =1;
    % leg = {};
    % names = {};
    % %custom legend
    % for i = ddvois
    %     leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.6);
    %     hold on
    %     names{c} = cst{i,2};
    %     c = c+1;
    % 
    % end
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.6);
    % names{c} = 'joint MinMax';
    % c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.6);
    % names{c} = 'joint SO';
    % c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.6);
    % names{c} = 'joint SU';
    % c = c+1;
    % % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
    % % names{c} = 'joint LET SU';
    % % c = c+1;
    % % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',2);
    % % names{c} = 'only Proton';
    % % c = c+1;
    % leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.6);
    % names{c} = 'Ref. Joint';
    % c = c+1;
    % axis([0,60,0,100])
    % legend ([leg{:}],names )
    % xlabel('Dirty BED [Gy]')
    % ax.YColor = 'none';
    % set(gca,'FontSize',18)


%% LETxDose DVH
DDtotalLETDVH  = matRad_calcDVH(cst_DDovernunder,result_DDovernunder{1,1}.mLETDose,'cum');
LETtotalLETDVH = matRad_calcDVH(cst_LETovernunder,result_LETovernunder{1,1}.mLETDose,'cum');
%onlyProtonDVH = matRad_calcDVH(cst,onlyProton{1}.mLETDose,'cum');
%refJointDVH = matRad_calcDVH(cst,RefJoint{1,1}.mLETDose,'cum');

subplot(1,2,2)
vois = [3];
for i  = vois
    plot(DDtotalLETDVH(i).doseGrid,DDtotalLETDVH(i).volumePoints,'LineWidth',1.2,'Color',cst_DDovernunder{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(LETtotalLETDVH(i).doseGrid,LETtotalLETDVH(i).volumePoints,'LineWidth',1.2,'Color',cst_LETovernunder{i,5}.visibleColor,'LineStyle', ':');
    %plot(onlyProtonDVH(i).doseGrid,onlyProtonDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    %plot(refJointDVH(i).doseGrid,refJointDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst_DDovernunder{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst_DDovernunder{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Joint DD';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'Joint LET';
c = c+1;
% leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.2);
% names{c} = 'only Proton';
% c = c+1;
% leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
% names{c} = 'Ref. Joint';
% c = c+1;
% axis([0,3,0,100])
legend ([leg{:}],names ) 
xlabel('LET')
ylabel('Volume [%]')
grid on


%% Quality indicators 

DDOU_physTotal= (pln_onlyPro(1).numOfFractions.*result_DDOU{1,1}.effect + pln_onlyPro(2).numOfFractions.*result_DDOU{1,2}.effect)./0.1;
DDOU_physProton = (pln_onlyPro(1).numOfFractions.*result_DDOU{1,1}.effect)/0.1;
DDOU_physPhoton = (pln_onlyPro(2).numOfFractions.*result_DDOU{1,2}.effect)/0.1;

LETOU_physTotal_2 = (pln_onlyPro(1).numOfFractions.*result_LETOU{1,1}.effect + pln_onlyPro(2).numOfFractions.*result_LETOU{1,2}.effect)./0.1;
LETOU_physProton = (pln_onlyPro(1).numOfFractions.*result_LETOU{1,1}.effect)/0.1;
LETOU_physPhoton = (pln_onlyPro(2).numOfFractions.*result_LETOU{1,2}.effect)/0.1;

DDOU_physQI =  matRad_calcQualityIndicators(cst,plnJO,DDOU_physTotal);
DDOU_proton_physQI =  matRad_calcQualityIndicators(cst,pln_onlyPro(1),DDOU_physProton);
DDOU_photon_physQI =  matRad_calcQualityIndicators(cst,pln_onlyPro(2),DDOU_physPhoton);


LETOU_physQI_2 =  matRad_calcQualityIndicators(cst,plnJO,LETOU_physTotal_2);
LETOU_proton_physQI = matRad_calcQualityIndicators(cst,pln_onlyPro(1),LETOU_physProton);
LETOU_photon_physQI = matRad_calcQualityIndicators(cst,pln_onlyPro(2),LETOU_physPhoton);


DDOU_physTotal= (pln_onlyPro(1).numOfFractions.*result_DDOU{1,1}.effect + pln_onlyPro(2).numOfFractions.*result_DDOU{1,2}.effect)./0.1;
LETOU_physTotal = (pln_onlyPro(1).numOfFractions.*result_LETOU{1,1}.effect + pln_onlyPro(2).numOfFractions.*result_LETOU{1,2}.effect)./0.1;

DDOU_physQI =  matRad_calcQualityIndicators(cst,plnJO,DDOU_physTotal);
LETOU_physQI =  matRad_calcQualityIndicators(cst,plnJO,LETOU_physTotal);

DDoverQI =  matRad_calcQualityIndicators(cst,plnJO,DDoverTotal)
LEToverQI =  matRad_calcQualityIndicators(cstQI,plnJO,LEToverTotal)
DDSOSUQI = matRad_calcQualityIndicators(cst,plnJO,DDOUDVH);
onlyProtonQI = matRad_calcQualityIndicators(cst,plnJO,onlyProtonTotal);
refJointQI = matRad_calcQualityIndicators(cst,plnJO,refJointTotal);

photonQI = matRad_calcQualityIndicators(cst,pln_onlyPro,photonPlan);

%% xtract Values
VOI = 7;
Quant = 'mean';

a(1) = onlyProtonQI(VOI).([Quant]);
a(2) = refJointQI(VOI).([Quant]);
a(3) = DDoverQI(VOI).([Quant]);
a(4) = LEToverQI(VOI).([Quant])

%% DD Quality indicators 
% LET 

DDoverQI =  matRad_calcQualityIndicators(cstQI,plnJO,result_DD_over{1,1}.dirtyDose)
LEToverQI =  matRad_calcQualityIndicators(cstQI,plnJO,result_LET_over{1,1}.dirtyDose)
onlyProtonQI = matRad_calcQualityIndicators(cst,plnJO, onlyProton{1}.dirtyDose)
refJointQI = matRad_calcQualityIndicators(cstQI,plnJO,RefJoint{1,1}.dirtyDose)
%% D BED 
n = 5;
a = 0.1;
b = 0.05 ;

%photon_DBED = effect( n,result_photon{1,1}.dirtyDose, a, b) ;

DDOU_DBED = effect( n,result_DDOU{1,1}.dirtyDose, a, b)/0.1 ;
DDOU_photonDBED = effect(25,result_DDOU{1,2}.dirtyDose,a,b)/0.1;

LETOU_DBED = effect( n,result_LETOU{1,1}.dirtyDose, a, b)/0.1 ;
LETOU_photonDBED = effect(25,result_LETOU{1,2}.dirtyDose,a,b)/0.1;

DDOU_BED = effect( n,result_DDOU{1,1}.physicalDose, a, b)/0.1 ;
LETOU_BED = effect( n,result_LETOU{1,1}.physicalDose, a, b)/0.1 ;

DDover_BED = effect( n,result_DD_over{1,1}.physicalDose, a, b)/0.1 ;

DDunder_DBED = effect( n,resultDD_Target{1,1}.dirtyDose, a, b)/0.1 ;
LETunder_DBED = effect( n,resultLET_Target{1,1}.dirtyDose, a, b)/0.1 ;
DDover_DBED = effect( n,result_DD_over{1,1}.dirtyDose, a, b) /0.1 ;
LETover_DBED = effect( n,result_LET_over{1,1}.dirtyDose, a, b)/0.1 ;
onlyProton_DBED = effect( 30,onlyProton{1}.dirtyDose, a, b)/0.1 ;
RefJoint_DBED = effect( n,RefJoint{1,1}.dirtyDose, a, b)/0.1 ;

%%
photon_DBED_QI = matRad_calcQualityIndicators(cst,pln_onlyPro,photon_DBED);

DDoverQI =  matRad_calcQualityIndicators(cstQI,plnJO,DD_DBED);
DDunder_DQI =  matRad_calcQualityIndicators(cst,plnJO,DDunder_DBED);

LETover_DQI =  matRad_calcQualityIndicators(cst,plnJO,LETover_DBED);
LETunder_DQI =  matRad_calcQualityIndicators(cst,plnJO,LETunder_DBED);

DDOU_DQI =  matRad_calcQualityIndicators(cst,plnJO,DDOU_DBED);
DDOU_photonDQI = matRad_calcQualityIndicators(cst,plnJO,DDOU_photonDBED);

LETOU_DQI =  matRad_calcQualityIndicators(cst,plnJO,LETOU_DBED);
LETOU_photonDQI = matRad_calcQualityIndicators(cst,plnJO,LETOU_photonDBED);

DDover_QI =  matRad_calcQualityIndicators(cst,plnJO,DDover_BED);


onlyProton_DQI = matRad_calcQualityIndicators(cst,plnJO,onlyProton_DBED );
refJoint_DQI = matRad_calcQualityIndicators(cst,plnJO,RefJoint_DBED);
%% 
VOI = 3;
Quant = 'std';

a(1) = onlyProtonQI(VOI).([Quant]);
a(2) = refJointQI(VOI).([Quant]);
a(3) = DDoverQI(VOI).([Quant]);
a(4) = LEToverQI(VOI).([Quant])


%% 
%% calc DVH for the current cube 

% load('FinalResultUnder.mat')
% load('D:\postDoc\LET_Paper\OnlyProton.mat')

DDunderTotal= (pln_onlyPro(1).numOfFractions.*resultDD_Target{1,1}.effect + pln_onlyPro(2).numOfFractions.*resultDD_Target{1,2}.effect)./0.1;
LETunderTotal = (pln_onlyPro(1).numOfFractions.*resultLET_Target{1,1}.effect + pln_onlyPro(2).numOfFractions.*resultLET_Target{1,2}.effect)./0.1;

onlyProtonTotal = effect(30,onlyProton{1,1}.RBExD,0.1,0.05)./0.1  ;
refJointTotal = (pln_onlyPro(1).numOfFractions.*RefJoint{1,1}.effect + pln_onlyPro(2).numOfFractions.*RefJoint{1,2}.effect)./0.1; 
% resultDVH  = matRad_calcDVH(cst,totalPlan,'cum');
% LET_DVH = resultDVH;
%% total BED DVH
cst = matRad_setOverlapPriorities(cst,ct.cubeDim);

DDtotalDVH  = matRad_calcDVH(cst,DDunderTotal,'cum');
LETtotalDVH = matRad_calcDVH(cst,LETunderTotal,'cum');
onlyProtonDVH = matRad_calcDVH(cst,onlyProtonTotal,'cum');
refJointDVH = matRad_calcDVH(cst,refJointTotal,'cum');

figure,
vois = [3,12,13,14,15];
for i  = vois
    plot(DDtotalDVH(i).doseGrid,DDtotalDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(LETtotalDVH(i).doseGrid,LETtotalDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    plot(onlyProtonDVH(i).doseGrid,onlyProtonDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    plot(refJointDVH(i).doseGrid,refJointDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Joint DD';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'Joint LET';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.2);
names{c} = 'only Proton';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
names{c} = 'Ref. Joint';
c = c+1;
xlabel('BED [Gy]')
ylabel('Volume [%]')
legend ([leg{:}],names ) 
set(gca,'FontSize',14)
grid on

%% Dirty Dose DVH
DDtotalDVH  = matRad_calcDVH(cst,resultDD_Target{1,1}.dirtyDose*5,'cum');
LETtotalDVH = matRad_calcDVH(cst,resultLET_Target{1,1}.dirtyDose*5,'cum');
onlyProtonDVH = matRad_calcDVH(cst,onlyProton{1}.dirtyDose*30,'cum');
refJointDVH = matRad_calcDVH(cst,RefJoint{1,1}.dirtyDose*5,'cum');


DDtotalDVH  = matRad_calcDVH(cst,resultDD_Target{1,1}.dirtyDose,'cum');
LETtotalDVH = matRad_calcDVH(cst,resultLET_Target{1,1}.dirtyDose,'cum');
onlyProtonDVH = matRad_calcDVH(cst,onlyProton{1}.dirtyDose,'cum');
refJointDVH = matRad_calcDVH(cst,RefJoint{1,1}.dirtyDose,'cum');

DDtotalDVH  = matRad_calcDVH(cst,effect(5,resultDD_Target{1,1}.dirtyDose,0.1,0.05),'cum');
LETtotalDVH = matRad_calcDVH(cst,effect(5,resultLET_Target{1,1}.dirtyDose ,0.1,0.05),'cum');
onlyProtonDVH = matRad_calcDVH(cst,effect(30, onlyProton{1}.dirtyDose,0.1,0.05),'cum');
refJointDVH = matRad_calcDVH(cst,effect(5,RefJoint{1,1}.dirtyDose ,0.1,0.05),'cum');
%%
%% D BED 
n = 5;
a = 0.1;
b = 0.05 ;

DD_DBED = effect( n,resultDD_Target{1,1}.dirtyDose, a, b) ;
LET_DBED = effect( n,resultLET_Target{1,1}.dirtyDose, a, b);
onlyProton_DBED = effect( 30,onlyProton{1}.dirtyDose, a, b);
RefJoint_DBED = effect( n,RefJoint{1,1}.dirtyDose, a, b);

photon_

DDoverQI =  matRad_calcQualityIndicators(cstQI,plnJO,DD_DBED);
LEToverQI =  matRad_calcQualityIndicators(cstQI,plnJO,LET_DBED);
onlyProtonQI = matRad_calcQualityIndicators(cstQI,plnJO,onlyProton_DBED );
refJointQI = matRad_calcQualityIndicators(cstQI,plnJO,RefJoint_DBED);


DDtotalDVH  = matRad_calcDVH(cst,DD_DBED,'cum');
LETtotalDVH = matRad_calcDVH(cst,LET_DBED,'cum');
onlyProtonDVH = matRad_calcDVH(cst,onlyProton_DBED,'cum');
refJointDVH = matRad_calcDVH(cst,RefJoint_DBED,'cum');

%%

fig = figure;

vois = [13];
for i  = vois
    plot(DDtotalDVH(i).doseGrid,DDtotalDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(LETtotalDVH(i).doseGrid,LETtotalDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    plot(onlyProtonDVH(i).doseGrid,onlyProtonDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    plot(refJointDVH(i).doseGrid,refJointDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '--');
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Joint DD';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'Joint LET';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',1.2);
names{c} = 'only Proton';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',1.2);
names{c} = 'Ref. Joint';
c = c+1;
% axis([0,3,0,100])
legend ([leg{:}],names ) 
xlabel(' Dirty BED [Gy]')
ylabel('Volume [%]')
set(gca,'FontSize',14)
grid on

%%
%% DD Quality indicators 
% LET 

DDoverQI =  matRad_calcQualityIndicators(cstQI,plnJO,resultDD_Target{1,1}.dirtyDose*5)
LEToverQI =  matRad_calcQualityIndicators(cstQI,plnJO,resultLET_Target{1,1}.dirtyDose*5)
onlyProtonQI = matRad_calcQualityIndicators(cstQI,plnJO, onlyProton{1}.dirtyDose*30)
refJointQI = matRad_calcQualityIndicators(cstQI,plnJO,RefJoint{1,1}.dirtyDose*5)

DDoverQI =  matRad_calcQualityIndicators(cstQI,plnJO,DDunderTotal)
LEToverQI =  matRad_calcQualityIndicators(cstQI,plnJO,LETunderTotal)
onlyProtonQI = matRad_calcQualityIndicators(cstQI,plnJO,onlyProtonTotal)
refJointQI = matRad_calcQualityIndicators(cstQI,plnJO,refJointTotal)

%% 
VOI = 13;
Quant = 'std';

a(1) = onlyProtonQI(VOI).([Quant]);
a(2) = refJointQI(VOI).([Quant]);
a(3) = DDoverQI(VOI).([Quant]);
a(4) = LEToverQI(VOI).([Quant])

%% QI Photon

PhotonQI0 =  matRad_calcQualityIndicators(cst,plnJO,result_photon{1,1}.physicalDose);

%%
a = PhotonQI(VOI).([Quant]);
