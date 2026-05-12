%% Master Thesis: box phantom comparisons based on --> carbon reference
matRad_rc;
% clear
matRad_cfg = MatRad_Config.instance();

load("BOXPHANTOM.mat")

%% set objectives
% use standard objectives

%% define models to compare

modelForComparison = ["LEMI","mMKM","aMCFMKM","MCFMKM"];

for i = 1:size(modelForComparison,2)
    cst_models_carbon.(modelForComparison(i)) = cst;
end

% deleting alphaX and betaX in phantom and setting model specific reference
% parameters

fn = fieldnames(cst_models_carbon);
bioParams = {'alphaX','betaX'};

for j = 1: numel(fn)
    for k = 1:size(cst_models_carbon.(fn{j}),1)
        if isfield(cst_models_carbon.(fn{j}){k,5},'alphaX') && isfield(cst_models_carbon.(fn{j}){k,5},'betaX')
            cst_models_carbon.(fn{j}){k,5} = rmfield(cst_models_carbon.(fn{j}){k,5},bioParams);
        end
        cst_models_carbon.(fn{j}){k,5}.bioParams = struct();
        cst_models_carbon.(fn{j}){k,5}.bioParams.cellLine = 'HSG';
        cst_models_carbon.(fn{j}){k,5}.bioParams.refRadiation = 'carbon';

        % take the alpha and beta parameters from the Inaniwa2015 box
        % phantom approach
        if strcmp(fn{j},'LEMI')
            cst_models_carbon.(fn{j}){k,5}.bioParams.alphaX = 0.9927;
            cst_models_carbon.(fn{j}){k,5}.bioParams.betaX = 0.0419;
        end

        % for models with different biological parameters
        if strcmp(fn{j},'mMKM')
            cst_models_carbon.(fn{j}){k,5}.bioParams.alphaX = 0.8878;
            cst_models_carbon.(fn{j}){k,5}.bioParams.betaX = 0.0615;
        end

        if strcmp(fn{j},'aMCFMKM')
            cst_models_carbon.(fn{j}){k,5}.bioParams.alphaX = 0.8733; 
            cst_models_carbon.(fn{j}){k,5}.bioParams.betaX = 0.0511;
        end
        if strcmp(fn{j},'MCFMKM')
            cst_models_carbon.(fn{j}){k,5}.bioParams.alphaX = 1.0018;
            cst_models_carbon.(fn{j}){k,5}.bioParams.betaX = 0.0500;
        end

    end
end


%% Pln Definition

for j = 1: numel(fn)

    pln_carbon.(fn{j}).radiationMode = 'carbon';
    %pln.(fn{j}).machine       = 'HITfluenceSpectra_updated';
    pln_carbon.(fn{j}).machine       = 'MCF_Fluence_updated_withoutRBE_origSingleNewSigma';
    
    %pln.(fn{j}).machine       = 'Generic_clusterDose_prestep_update';
    pln_carbon.(fn{j}).propDoseCalc.calcLET = 0;

    pln_carbon.(fn{j}).numOfFractions       = 30;
    pln_carbon.(fn{j}).propStf.gantryAngles = 90; % for the box phantom
    pln_carbon.(fn{j}).propStf.couchAngles   = zeros(numel(pln_carbon.(fn{j}).propStf.gantryAngles),1);
    pln_carbon.(fn{j}).propStf.bixelWidth    = 3;
    pln_carbon.(fn{j}).propStf.numOfBeams    = numel(pln_carbon.(fn{j}).propStf.gantryAngles);

    pln_carbon.(fn{j}).propStf.isoCenter     = ones(pln_carbon.(fn{j}).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_models_carbon.(fn{j}),ct,0);
    pln_carbon.(fn{j}).propStf.longitudinalSpotSpacing = 1;
    pln_carbon.(fn{j}).propOpt.runDAO        = 0;
    pln_carbon.(fn{j}).propSeq.runSequencing = 0;

    % dose calculation settings
    pln_carbon.(fn{j}).propDoseCalc.doseGrid.resolution.x = 3;
    pln_carbon.(fn{j}).propDoseCalc.doseGrid.resolution.y = 3;
    pln_carbon.(fn{j}).propDoseCalc.doseGrid.resolution.z = 3;

    pln_carbon.(fn{j}).multScen = matRad_multScen(ct, 'nomScen');
end
%%
% pln_2.mMKM_C = pln_2.mMKM;
% pln_2.MCFMKM_C = pln_2.MCFMKM;
% pln_2.MCFMKM_C.machine = 'MCF_Fluence_updated_double_carbon';

%% stf generation
for j = 1: numel(fn)
    stf_carbon.(fn{j}) = matRad_generateStf(ct,cst_models_carbon.(fn{j}),pln_carbon.(fn{j}));
end
%stf = matRad_generateSingleBixelStf(ct,cst,pln);

%% Pln RBE model

% change the names of the RBEtables before running it!

for j = 1: numel(fn)
    pln_carbon.(fn{j}).bioModel = matRad_bioModel(pln_carbon.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');

    if strcmp(fn{j},'LEMI')
        pln_carbon.(fn{j}).bioModel.quantityTableName      = 'RBEtable_LEMI_carbon_217_MCF_AX9927_BX419';

    end

    if strcmp(fn{j},'mMKM')
        pln_carbon.(fn{j}).bioModel.quantityTableName = 'RBEtable_mMKM_carbon_A0_217_AX8878_BX615';

    end

    if strcmp(fn{j},'aMCFMKM')
        pln_carbon.(fn{j}).bioModel.quantityTableName = 'RBEtable_approxMCFMKM_carbon_AX8733_BX511';

    end

    if strcmp(fn{j},'MCFMKM')
        pln_carbon.(fn{j}).bioModel.quantityTableName = 'RBEtable_MCFMKM_Shannon_carbon_AX10018_BX005';

    end


    pln_carbon.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

    selectZ = [1,2,3,4,5,6];

    selectA = [1,4,7,9,11,12];

    pln_carbon.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);

    pln_carbon.(fn{j}).propOpt.quantityOpt = 'physicalDose';
    pln_carbon.(fn{j}).propDoseCalc.engine = 'HongPB';
end

%% Dose calculation

% everything is based on the LEMI optimization (weight vector)
% only calculate the dij once for LEMI

dij_carbon.LEMI = matRad_calcDoseInfluence(ct,cst_models_carbon.(fn{1}),stf.(fn{1}),pln_carbon.(fn{1}));

%% Optimization

% set optimization quantity
% for the box phantom the quantity should be physical dose, so the dose
% stays the same and the RBE can be analyzed

for j = 1: numel(fn)
    for k = 1:size(cst_models_carbon.(fn{j}),1)
        if isempty(cst_models_carbon.(fn{j}){k,6})
            cst_models_carbon.(fn{j}){k,6} = {};
        else
            cst_models_carbon.(fn{j}){k,6}{1}.quantity = 'physicalDose';
        end
    end
end

resultGUI_carbon.LEMI = matRad_fluenceOptimization(dij_carbon.LEMI,cst_models_carbon.LEMI,pln_carbon.LEMI);

%% Recalculation
% same weights for the same dose for every biological model

for j = 2:numel(fn)
    resultGUI_carbon.(fn{j}) = matRad_calcDoseForward(ct,cst_models_carbon.(fn{j}),stf_carbon.(fn{j}),pln_carbon.(fn{j}),resultGUI_carbon.LEMI.w);
end

%% Plot the SOBP

f = figure('WindowState','maximized');
subplot(2,3,1)

blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
purple = [1 0.1840 0.5560];
brown = [0.6350 0.0780 0.1840];
blueD = [0 0 1];
colors = [blue; orange; green; red; purple; brown; blueD];

lineWidth = 2;
markerSize = 15;
yIndex   = 80;
x = dij.LEMI.ctGrid.x;
x = x(41:117);
% look at the plot physicalDose in depth before cutting the begin and the
% ending
% x(x < -120) = [];
% x(x > 108) = [];
x = x + 120;

for j = 1:numel(fn)
    profile_PB_d.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).physicalDose(yIndex,:,yIndex))));
    plot(x,profile_PB_d.(fn{j})(41:117),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('physicalDose [Gy]','FontSize',20);
axis([0,230,0,2.5])
ax = gca();

legend('FontSize',20);

subplot(2,3,2)

for j = 1:numel(fn)
    profile_PB_r.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_PB_r.(fn{j})(41:117),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('RBExDose [Gy(RBE)]','FontSize',20);
axis([0,230,0,3.5])
ax = gca();

legend('FontSize',20);

subplot(2,3,3)

for j = 1:numel(fn)
    profile_PB_e.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).effect(yIndex,:,yIndex))));
    plot(x,profile_PB_e.(fn{j})(41:117),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('effect','FontSize',20);
axis([0,230,0,3.5])
ax = gca();

legend('FontSize',20);

subplot(2,3,4)

for j = 1:numel(fn)
    profile_PB_a.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).alphaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_a.(fn{j})(41:117),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('alphaDose','FontSize',20);
axis([0,230,0,3.5])
ax = gca();

legend('FontSize',20);

subplot(2,3,5)

for j = 1:numel(fn)
    profile_PB_b.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).SqrtBetaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_b.(fn{j})(41:117),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('SqrtBetaDose','FontSize',20);
axis([0,230,0,0.6])
ax = gca();

legend('FontSize',20);


subplot(2,3,6)

for j = 1:numel(fn)
    BED.(fn{j}).BED = resultGUI.(fn{j}).effect ./ cst_models_carbon.(fn{j}){k,5}.bioParams.alphaX;
    profile_PB_bed.(fn{j}) = flip((squeeze(BED.(fn{j}).BED(yIndex,:,yIndex))));
    plot(x,profile_PB_bed.(fn{j})(41:117),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('BED [Gy]','FontSize',20);
axis([0,230,0,4])
ax = gca();

legend('FontSize',20);
