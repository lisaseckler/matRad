%% Master Thesis: TG166 prostate comparisons based on --> photon reference
%matRad_rc;
clear
matRad_cfg = MatRad_Config.instance();

load("PROSTATE.mat")

%% set objectives
cst{9,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,20));
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));
cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cst{6,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(1000,52));
cst{7,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(1000,40));

%% define models to compare

modelForComparison = ["LEMI","mMKM","aMCFMKM","MCFMKM"];

for i = 1:size(modelForComparison,2)
    cst_models_photon.(modelForComparison(i)) = cst;
end

% deleting alphaX and betaX in phantom and setting model specific reference
% parameters

fn = fieldnames(cst_models_photon);
bioParams = {'alphaX','betaX'};

for j = 1: numel(fn)
    for k = 1:size(cst_models_photon.(fn{j}),1)
        if isfield(cst_models_photon.(fn{j}){k,5},'alphaX') && isfield(cst_models_photon.(fn{j}){k,5},'betaX')
            cst_models_photon.(fn{j}){k,5} = rmfield(cst_models_photon.(fn{j}){k,5},bioParams);
        end
        cst_models_photon.(fn{j}){k,5}.bioParams = struct();
        cst_models_photon.(fn{j}){k,5}.bioParams.cellLine = 'HSG';
        cst_models_photon.(fn{j}){k,5}.bioParams.refRadiation = 'photons'; 
        cst_models_photon.(fn{j}){k,5}.bioParams.alphaX = 0.217;
        cst_models_photon.(fn{j}){k,5}.bioParams.betaX = 0.0615;
    end
end


%% Pln Definition

for j = 1: numel(fn)

    pln_photon.(fn{j}).radiationMode = 'carbon';
    pln_photon.(fn{j}).machine       = 'HITfluenceSpectra_updatedBig';
    pln_photon.(fn{j}).propDoseCalc.calcLET = 0;

    pln_photon.(fn{j}).numOfFractions       = 12;
    pln_photon.(fn{j}).propStf.gantryAngles = [90 270]; %[45 135 225 315];
    pln_photon.(fn{j}).propStf.couchAngles   = zeros(numel(pln_photon.(fn{j}).propStf.gantryAngles),1);
    pln_photon.(fn{j}).propStf.bixelWidth    = 5;
    pln_photon.(fn{j}).propStf.numOfBeams    = numel(pln_photon.(fn{j}).propStf.gantryAngles);

    pln_photon.(fn{j}).propStf.isoCenter     = ones(pln_photon.(fn{j}).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_models_photon.(fn{j}),ct,0);
    pln_photon.(fn{j}).propStf.longitudinalSpotSpacing = 1;
    pln_photon.(fn{j}).propOpt.runDAO        = 0;
    pln_photon.(fn{j}).propSeq.runSequencing = 0;

    % dose calculation settings
    pln_photon.(fn{j}).propDoseCalc.doseGrid.resolution.x = 3;
    pln_photon.(fn{j}).propDoseCalc.doseGrid.resolution.y = 3;
    pln_photon.(fn{j}).propDoseCalc.doseGrid.resolution.z = 3;


    pln_photon.(fn{j}).multScen = matRad_multScen(ct, 'nomScen');
end

%% stf generation
for j = 1:numel(fn)
    stf_photon.(fn{j}) = matRad_generateStf(ct,cst_models_photon.(fn{j}),pln_photon.(fn{j}));
end
%stf.LEMI = matRad_generateSingleBixelStf(ct,cst_models.LEMI,pln.LEMI);

%% Pln RBE model

for j = 1: numel(fn)

    if strcmp(fn{j},'LEMI')
        pln_photon.(fn{j}).bioModel = matRad_bioModel(pln_photon.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln_photon.(fn{j}).bioModel.quantityTableName      = 'RBEtable_LEMI_AX217_BX615';
        pln_photon.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln_photon.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);


    end

    if strcmp(fn{j},'mMKM')
        pln_photon.(fn{j}).bioModel = matRad_bioModel(pln_photon.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');

        pln_photon.(fn{j}).bioModel.quantityTableName = 'RBEtable_mMKM_photon_A0_AX217_BX615';
        pln_photon.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln_photon.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);


    end
    if strcmp(fn{j},'aMCFMKM')
        pln_photon.(fn{j}).bioModel = matRad_bioModel(pln_photon.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln_photon.(fn{j}).bioModel.quantityTableName = 'RBEtable_MCFMKM_photon_AX0217_BX00615';
        pln_photon.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln_photon.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);

    end
    if strcmp(fn{j},'MCFMKM')
        pln_photon.(fn{j}).bioModel = matRad_bioModel(pln_photon.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln_photon.(fn{j}).bioModel.quantityTableName = 'RBEtable_MCFMKM_photon_Shannon_AX117_BX615';
        pln_photon.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln_photon.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);

    end



    pln_photon.(fn{j}).propOpt.quantityOpt = 'physicalDose';
    pln_photon.(fn{j}).propDoseCalc.engine = 'HongPB';
end

%% Dose calculation

% everything is based on the LEMI optimization (weight vector)
% only calculate the dij once for LEMI

% TabulatedDoseAveragedKernelModel is modified: dE/dx is transposed because
% of MCF_Fluence machine --> denumerator and numerator
dij_photon.LEMI = matRad_calcDoseInfluence(ct,cst_models_photon.(fn{1}),stf_photon.(fn{1}),pln_photon.(fn{1}));

%% Optimization

% set optimization quantity
% for the box phantom the quantity should be physical dose, so the dose
% stays the same and the RBE can be analyzed

for j = 1: numel(fn)
    for k = 1:size(cst_models_photon.(fn{j}),1)
        if isempty(cst_models_photon.(fn{j}){k,6})
            cst_models_photon.(fn{j}){k,6} = {};
        else
            cst_models_photon.(fn{j}){k,6}{1}.quantity = 'effect';
        end
    end
end

resultGUI_testPhoton.LEMI = matRad_fluenceOptimization(dij_test.LEMI,cst_models_photon.LEMI,pln_photon.LEMI);

%% TOPAS
pln_TOPASphoton_Test = pln_photon.LEMI;
%numel(fn)
pln_TOPASphoton_Test.propDoseCalc.engine = 'TOPAS';

% if pln_TOPAS.propDoseCalc.engine == 'TOPAS'
%     if pln.bioModel.model == 'doseAveragedTabulatedAlphaBeta'
%         pln_TOPAS = pln;
%         pln_TOPAS.propDoseCalc.engine = 'TOPAS';
%         pln_TOPAS.bioModel = matRad_bioModel(pln.radiationMode, 'doseAveragedTabulatedAlphaBeta');
%         pln_TOPAS.bioModel.RBEtableName = 'RBEtable_MCFMKM_0313';
%         pln_TOPAS.bioModel.includedFragments = struct('Z',[1,2,3,4,5,6], ...
%                                         'A',[1,4,7,9,11,12]);
%         %pln_TOPAS.bioModel = pln_TOPAS.bioModel.setFragmentIndexesInBaseData(pln.bioModel.fragmentIndexesInBaseData);
%         pln_TOPAS.bioModel = pln_TOPAS.bioModel.setTissueAlphaX(pln.bioModel.quantityTable.meta.alphaX);
%         pln_TOPAS.bioModel = pln_TOPAS.bioModel.setTissueBetaX(pln.bioModel.quantityTable.meta.betaX);
%
%     end
% end


pln_TOPASphoton_Test.propDoseCalc.externalCalculation = 'write';
pln_TOPASphoton_Test.propDoseCalc.calcBioDose = true;
pln_TOPASphoton_Test.propDoseCalc.numHistoriesDirect = 1000000;
pln_TOPASphoton_Test.propDoseCalc.numOfRuns = 5;

%resultGUI_TOPAS.(fn{j}) = matRad_calcDoseForward(ct,cst_models.(fn{j}),stf_carbon.(fn{j}),pln_TOPASphoton.(fn{j}),resultGUI_effect.LEMI.w);


%%
resultGUI_TOPASphoton_Test.LEMI = matRad_calcDoseForward(ct,cst_models_photon.LEMI,stf_test.LEMI,pln_TOPASphoton_Test,resultGUI_testPhoton.LEMI.w);

%% if TOPAS calculation is done
pln_TOPASphoton_Test.bioModel = pln_carbon.mMKM.bioModel;
pln_TOPASphoton_Test.propDoseCalc.externalCalculation = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\userdata\TOPAS\Results\carbon_HITfluenceSpectra_updatedBig_prostatemMKMcarbon_19-03-26';
resultGUI_fullTOPAS_carbon.mMKM = matRad_calcDoseForward(ct,cst_models_carbon.mMKM,stf_carbon.mMKM,pln_TOPASphoton_Test,resultGUI_carbon.LEMI.w);

%% delete dijs so far
clearvars -regexp ^dij


%% Recalculation
% same weights for the same dose for every biological model

for j = 2:numel(fn)
    resultGUI.(fn{j}) = matRad_calcDoseForward(ct,cst_models_photon.(fn{j}),stf.(fn{j}),pln_photon.(fn{j}),resultGUI.LEMI.w);
end

%% plot Slice Wrapper

plane = 3;
slice = 25;

LEMI_physD = resultGUI_fullTOPAS.LEMI.physicalDose;
% LEMII_physD = resultGUI_effect.LEMII.physicalDose;
% LEMIII_physD = resultGUI_effect.LEMIII.physicalDose;
% MKM_physD = resultGUI_effect.MKM.physicalDose;
mMKM_physD = resultGUI_fullTOPAS.mMKM.physicalDose;
aMCFMKM_physD = resultGUI_fullTOPAS.aMCFMKM.physicalDose;
MCFMKM_physD = resultGUI_fullTOPAS.MCFMKM.physicalDose;

LEMI_RBExD = resultGUI_fullTOPAS.LEMI.RBExDose;
% LEMII_RBExD = resultGUI_effect.LEMII.RBExDose;
% LEMIII_RBExD = resultGUI_effect.LEMIII.RBExDose;
% MKM_RBExD = resultGUI_effect.MKM.RBExDose;
mMKM_RBExD = resultGUI_fullTOPAS.mMKM.RBExDose;
aMCFMKM_RBExD = resultGUI_fullTOPAS.aMCFMKM.RBExDose;
MCFMKM_RBExD = resultGUI_fullTOPAS.MCFMKM.RBExDose;

LEMIc_RBExD = resultGUI_fullTOPAS_carbon.LEMI.RBExDose ./ resultGUI_fullTOPAS_carbon.LEMI.physicalDose .* resultGUI_fullTOPAS.LEMI.physicalDose;
% LEMII_RBExD = resultGUI_effect.LEMII.RBExDose;
% LEMIII_RBExD = resultGUI_effect.LEMIII.RBExDose;
% MKM_RBExD = resultGUI_effect.MKM.RBExDose;
mMKMc_RBExD = resultGUI_fullTOPAS_carbon.mMKM.RBExDose ./ resultGUI_fullTOPAS_carbon.mMKM.physicalDose .* resultGUI_fullTOPAS.mMKM.physicalDose;
aMCFMKMc_RBExD = resultGUI_fullTOPAS_carbon.aMCFMKM.RBExDose ./ resultGUI_fullTOPAS_carbon.aMCFMKM.physicalDose .* resultGUI_fullTOPAS.aMCFMKM.physicalDose;
MCFMKMc_RBExD = resultGUI_fullTOPAS_carbon.MCFMKM.RBExDose ./ resultGUI_fullTOPAS_carbon.MCFMKM.physicalDose .* resultGUI_fullTOPAS.MCFMKM.physicalDose;

LEMI_effect = resultGUI_fullTOPAS.LEMI.effect;
% LEMII_effect = resultGUI_effect.LEMII.effect;
% LEMIII_effect = resultGUI_effect.LEMIII.effect;
% MKM_effect = resultGUI_effect.MKM.effect;
mMKM_effect = resultGUI_fullTOPAS.mMKM.effect;
aMCFMKM_effect = resultGUI_fullTOPAS.aMCFMKM.effect;
MCFMKM_effect = resultGUI_fullTOPAS.MCFMKM.effect;

%%
f = figure('WindowState','maximized');

t = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');

set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

set(groot,'defaultAxesFontSize',14)
set(groot,'defaultTextFontSize',14)

%subplot(3,4,1)
nexttile;
doseWindow = [0 2];

matRad_plotSliceWrapper(gca,ct,cst,1,LEMI_physD,plane,slice,[],[],colorcube,[],doseWindow);
title('LEM-I')
xlabel('')
xticklabels([])
yticklabels([])
colorbar off
text(0.95,0.95,'(a)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
ylabel('$D_{phys}$ [Gy]')
zoom(2)

%subplot(3,4,2)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,mMKM_physD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
colorbar off
text(0.95,0.95,'(b)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
title('mMKM')
zoom(2)

%subplot(3,4,3)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,aMCFMKM_physD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
colorbar off
text(0.95,0.95,'(c)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
title('approx. MCF-MKM')
zoom(2)

%subplot(3,4,4)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,MCFMKM_physD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
text(0.95,0.95,'(d)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
title('MCF-MKM')
zoom(2)

%subplot(3,4,5)
nexttile;
doseWindow = [0 5];

matRad_plotSliceWrapper(gca,ct,cst,1,LEMI_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
title('')
colorbar off
text(0.95,0.95,'(e)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
ylabel('$D_{bio,\gamma}$ [Gy(RBE)]')
zoom(2)

%subplot(3,4,6)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,mMKM_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
colorbar off
text(0.95,0.95,'(f)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)

%subplot(3,4,7)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,aMCFMKM_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
colorbar off
text(0.95,0.95,'(g)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)

%subplot(3,4,8)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,MCFMKM_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
text(0.95,0.95,'(h)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)

nexttile;
doseWindow = [0 2.5];

matRad_plotSliceWrapper(gca,ct,cst,1,LEMIc_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
title('')
colorbar off
text(0.95,0.95,'(i)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
ylabel('$D_{bio,c}$ [Gy(RBE)]')
zoom(2)

%subplot(3,4,6)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,mMKMc_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
colorbar off
text(0.95,0.95,'(j)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)

%subplot(3,4,7)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,aMCFMKMc_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
colorbar off
text(0.95,0.95,'(k)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)

%subplot(3,4,8)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,MCFMKMc_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
text(0.95,0.95,'(l)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)


%subplot(3,4,9)
nexttile;
doseWindow = [0 3];

matRad_plotSliceWrapper(gca,ct,cst,1,LEMI_effect,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
title('')
colorbar off
text(0.95,0.95,'(m)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
ylabel('effect')
zoom(2)

%subplot(3,4,10)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,mMKM_effect,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
colorbar off
text(0.95,0.95,'(n)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)

%subplot(3,4,11)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,aMCFMKM_effect,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
colorbar off
text(0.95,0.95,'(o)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)

%subplot(3,4,12)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,MCFMKM_effect,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
text(0.95,0.95,'(p)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)