% box phantom comparisons based on --> photon reference
%matRad_rc;
clear
matRad_cfg = MatRad_Config.instance();

load("BOXPHANTOM.mat")

%% set objectives
% use standard objectives
% penalty from tumor set to 1000

%% define models to compare

modelForComparison = ["LEMI","mMKM","MCFMKM"];
cst{2,6}{1,1}.penalty = 1500;

for i = 1:size(modelForComparison,2)
    cst_models.(modelForComparison(i)) = cst;
end

% deleting alphaX and betaX in phantom and setting model specific reference
% parameters

fn = fieldnames(cst_models);
bioParams = {'alphaX','betaX'};

for j = 1:numel(fn)
    for k = 1:size(cst_models.(fn{j}),1)
        if isfield(cst_models.(fn{j}){k,5},'alphaX') && isfield(cst_models.(fn{j}){k,5},'betaX')
            cst_models.(fn{j}){k,5} = rmfield(cst_models.(fn{j}){k,5},bioParams);
        end
        cst_models.(fn{j}){k,5}.bioParams = struct();
        cst_models.(fn{j}){k,5}.bioParams.cellLine = 'HSG';
        cst_models.(fn{j}){k,5}.bioParams.refRadiation = 'photons';
        if strcmp(fn{j},'LEMI')
            cst_models.(fn{j}){k,5}.bioParams.alphaX = 0.313;
            cst_models.(fn{j}){k,5}.bioParams.betaX = 0.0615;
        end
        if strcmp(fn{j},'mMKM')
            cst_models.(fn{j}){k,5}.bioParams.alphaX = 0.3312;
            cst_models.(fn{j}){k,5}.bioParams.betaX = 0.0615;

        end
        if strcmp(fn{j},'MCFMKM')
            cst_models.(fn{j}){k,5}.bioParams.alphaX = 0.217;
            cst_models.(fn{j}){k,5}.bioParams.betaX = 0.0615;

        end
        % cst_models.(fn{j}){k,5}.bioParams.refRadiation = 'carbon';
        % if strcmp(fn{j},'LEMI')
        %     cst_models.(fn{j}){k,5}.bioParams.alphaX = 0.97;
        %     cst_models.(fn{j}){k,5}.bioParams.betaX = 0.0419;
        % end
        % if strcmp(fn{j},'mMKM')
        %     cst_models.(fn{j}){k,5}.bioParams.alphaX = 0.7554;
        %     cst_models.(fn{j}){k,5}.bioParams.betaX = 0.0615;
        % 
        % end
        % if strcmp(fn{j},'MCFMKM')
        %     cst_models.(fn{j}){k,5}.bioParams.alphaX = 1.019;
        %     cst_models.(fn{j}){k,5}.bioParams.betaX = 0.0498;
        % 
        % end
    end
end


%% Pln Definition

for j = 1:numel(fn)% numel(fn)

    pln.(fn{j}).radiationMode = 'carbon';
    %pln.(fn{j}).machine       = 'MCF_doubleGauss_ridgeFilter_sigmaAveraged';
    pln.(fn{j}).machine        = 'HITfluenceSpectra_withEnergySpectrum';
    pln.(fn{j}).propDoseCalc.calcLET = 0;

    pln.(fn{j}).numOfFractions       = 30;
    pln.(fn{j}).propStf.gantryAngles = 90; % for the box phantom
    pln.(fn{j}).propStf.couchAngles   = zeros(numel(pln.(fn{j}).propStf.gantryAngles),1);
    pln.(fn{j}).propStf.bixelWidth    = 5;
    pln.(fn{j}).propStf.numOfBeams    = numel(pln.(fn{j}).propStf.gantryAngles);

    pln.(fn{j}).propStf.isoCenter     = ones(pln.(fn{j}).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_models.(fn{j}),ct,0);
    pln.(fn{j}).propStf.longitudinalSpotSpacing = 1;
    pln.(fn{j}).propOpt.runDAO        = 0;
    pln.(fn{j}).propSeq.runSequencing = 0;

    % dose calculation settings
    pln.(fn{j}).propDoseCalc.doseGrid.resolution.x = 3;
    pln.(fn{j}).propDoseCalc.doseGrid.resolution.y = 3;
    pln.(fn{j}).propDoseCalc.doseGrid.resolution.z = 3;

    pln.(fn{j}).multScen = matRad_multScen(ct, 'nomScen');
end

%% stf generation
for j = 1: numel(fn)
    stf.(fn{j}) = matRad_generateStf(ct,cst_models.(fn{j}),pln.(fn{j}));
end

%% Pln RBE model

for j = 1:numel(fn)% numel(fn)

    if strcmp(fn{j},'LEMI')
        pln.(fn{j}).bioModel = matRad_bioModel(pln.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln.(fn{j}).bioModel.quantityTableName      = 'RBEtable_LEMI_Scholz06_AX313_BX615';

        pln.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);


    end

    if strcmp(fn{j},'mMKM')
        pln.(fn{j}).bioModel = matRad_bioModel(pln.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');

        pln.(fn{j}).bioModel.quantityTableName = 'RBEtable_mMKM_Inaniwa10_A0B0_AX3312_BX593';
        
        pln.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);

    end

    if strcmp(fn{j},'MCFMKM')
        pln.(fn{j}).bioModel = matRad_bioModel(pln.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln.(fn{j}).bioModel.quantityTableName = 'RBEtable_MCFMKM_photon_Shannon_AX117_BX615';
        pln.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);

    end


    pln.(fn{j}).bioModel.calcLET = 1;

    pln.(fn{j}).propOpt.quantityOpt = 'physicalDose';
    pln.(fn{j}).propDoseCalc.engine = 'HongPB';
end

%% Dose calculation

dij.MCFMKM = matRad_calcDoseInfluence(ct,cst_models.(fn{3}),stf.(fn{3}),pln.(fn{3}));

%% Optimization

% set optimization quantity
% for the box phantom the quantity should be physical dose, so the dose
% stays the same and the RBE can be analyzed

for j = 1: numel(fn)
    for k = 1:size(cst_models.(fn{j}),1)
        if isempty(cst_models.(fn{j}){k,6})
            cst_models.(fn{j}){k,6} = {};
        else
            cst_models.(fn{j}){k,6}{1}.quantity = 'physicalDose';
        end
    end
end

resultGUI.MCFMKM = matRad_fluenceOptimization(dij.(fn{3}),cst_models.(fn{3}),pln.(fn{3}));

%% Recalculation
% same weights for the same dose for every biological model

for j = [1]
    resultGUI_carbon_MCF.(fn{j}) = matRad_calcDoseForward(ct,cst_models_carbon.(fn{j}),stf_MCF.(fn{j}),pln_MCF.(fn{j}),resultGUI_carbon_MCF.(fn{3}).w);
end

%% Plot the SOBP

f = figure('WindowState','maximized');

t = tiledlayout(2,3,'TileSpacing','tight','Padding','loose');

set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

set(groot,'defaultAxesFontSize',14)
set(groot,'defaultTextFontSize',14)

%subplot(2,3,1)

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
x = dij_HIT.MCFMKM.ctGrid.x;
% look at the plot physicalDose in depth before cutting the begin and the
% ending
x = x(41:121);
% x(x < -120) = [];
% x(x > 123) = [];
x = x + 120;

nexttile;
for j = [1 2 3]
    profile_PB_d.(fn{j}) = flip((squeeze(resultGUI_carbon_HIT.(fn{j}).physicalDose(yIndex,:,yIndex))));
    plot(x,profile_PB_d.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_d.(fn{j}) = flip((squeeze(resultGUI_carbon_MCF.(fn{j}).physicalDose(yIndex,:,yIndex))));
    plot(x,profile_MCF_d.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
% xlabel('depth [mm]', 'FontSize',20);
xlabel('')
xticklabels([])
ylabel('$D_{phys}$ [Gy]');
text(0.95,0.95,'(a)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,2.5])
ax = gca();

%subplot(2,3,2)
nexttile;

for j = [1 2 3]
    profile_PB_r.(fn{j}) = flip((squeeze(resultGUI_carbon_HIT.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_PB_r.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_r.(fn{j}) = flip((squeeze(resultGUI_carbon_MCF.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_MCF_r.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
%xlabel('depth [mm]', 'FontSize',20);
xlabel('')
xticklabels([])
ylabel('$D_{bio,c}$ [Gy(RBE)]');
text(0.95,0.95,'(b)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,3])
ax = gca();

% nexttile;
% 
% for j = [1 2 3]
%     profile_PB_rBc.(fn{j}) = flip((squeeze(resultGUI_carbon_HIT.(fn{j}).RBExDose(yIndex,:,yIndex))));
%     plot(x,profile_PB_rBc.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
%     hold on
%     profile_MCF_rBc.(fn{j}) = flip((squeeze(resultGUI_carbon_MCF.(fn{j}).RBExDose(yIndex,:,yIndex))));
%     plot(x,profile_MCF_rBc.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
%     hold on
% end
% 
% 
% grid on;
% xlabel('')
% xticklabels([])
% ylabel('$D_{bio,c}$ [Gy(RBE)]');
% text(0.95,0.95,'(c)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
% axis([0,250,0,3.5])
% ax = gca();

%subplot(2,3,4)
nexttile;
for j = [1 2 3]
    profile_PB_a.(fn{j}) = flip((squeeze(resultGUI_carbon_HIT.(fn{j}).alphaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_a.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_a.(fn{j}) = flip((squeeze(resultGUI_carbon_MCF.(fn{j}).alphaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_MCF_a.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('depth [mm]');
ylabel('$\alpha \times$ dose');
text(0.95,0.95,'(c)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,3.5])
ax = gca();

%subplot(2,3,5)
nexttile;
for j = [1 2 3]
    profile_PB_b.(fn{j}) = flip((squeeze(resultGUI_carbon_HIT.(fn{j}).SqrtBetaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_b.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_b.(fn{j}) = flip((squeeze(resultGUI_carbon_MCF.(fn{j}).SqrtBetaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_MCF_b.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('depth [mm]');
ylabel('$\sqrt{\beta} \times$ dose');
text(0.95,0.95,'(d)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,0.6])
ax = gca();

%subplot(2,3,6)
nexttile;

for j = [1 2 3]
    profile_PB_e.(fn{j}) = flip((squeeze(resultGUI_carbon_HIT.(fn{j}).effect(yIndex,:,yIndex))));
    plot(x,profile_PB_e.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_e.(fn{j}) = flip((squeeze(resultGUI_carbon_MCF.(fn{j}).effect(yIndex,:,yIndex))));
    plot(x,profile_MCF_e.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('depth [mm]');
ylabel('effect');
text(0.95,0.95,'(e)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,4])
ax = gca();

% Legend
hold on
% --- Dummy handles for MODEL (color) ---
cLEM = plot(nan, nan, '-', 'Color', colors(1,:), 'LineWidth', lineWidth); % j=1, blue
cMKM = plot(nan, nan, '-', 'Color', colors(2,:), 'LineWidth', lineWidth); % j=2, orange
cMCF = plot(nan, nan, '-', 'Color', colors(3,:), 'LineWidth', lineWidth); % j=3, green

% --- Dummy handles for MACHINE TYPE (linestyle) ---
sHIT = plot(nan, nan, '-',  'Color', 'k', 'LineWidth', lineWidth); % solid = HIT
sMCF = plot(nan, nan, '--', 'Color', 'k', 'LineWidth', lineWidth); % dashed = MCF

lgd = legend([cLEM cMKM cMCF sHIT sMCF], ...
    {'LEM-I','mMKM','MCF-MKM','HIT','MCF'}, ...
    'Interpreter','latex', ...
    'Orientation','horizontal', ...
    'NumColumns',3);

lgd.Layout.Tile = 'south';

set(gcf,'Position',[100 100 900 500])

