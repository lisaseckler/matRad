%% Box phantom comparisons between HIT and MCF based on --> photon reference
%matRad_rc;
clear
matRad_cfg = MatRad_Config.instance();

load("BOXPHANTOM.mat")

%% define models to compare

modelForComparison = ["LEMI"];

% increase tumor penalty for better dose homogenity
cst{2,6}{1,1}.penalty = 1500;

for i = 1:size(modelForComparison,2)
    cst_models.(modelForComparison(i)) = cst;
end

% deleting alphaX and betaX in phantom and setting model specific reference
% parameters

fn = fieldnames(cst_models);
bioParams = {'alphaX','betaX'};

% Photon Reference

for j = 1: numel(fn)
    for k = 1:size(cst_models.(fn{j}),1)
        if isfield(cst_models.(fn{j}){k,5},'alphaX') && isfield(cst_models.(fn{j}){k,5},'betaX')
            cst_models.(fn{j}){k,5} = rmfield(cst_models.(fn{j}){k,5},bioParams);
        end
        cst_models.(fn{j}){k,5}.bioParams = struct();
        cst_models.(fn{j}){k,5}.bioParams.cellLine = 'HSG';
        cst_models.(fn{j}){k,5}.bioParams.refRadiation = 'photons'; 
        cst_models.(fn{j}){k,5}.bioParams.alphaX = 0.217;
        cst_models.(fn{j}){k,5}.bioParams.betaX = 0.0615;

    end
end

%% Pln Definition

% change the index between MCF and HIT to make sure not to overwrite your
% variables

for j = 1:numel(fn)

    pln_HIT.(fn{j}).radiationMode = 'carbon';
    pln_HIT.(fn{j}).machine        = 'HITfluenceSpectra_updatedBig_singleGaussianFit_v2';

    pln_HIT.(fn{j}).propDoseCalc.calcLET = 0; % for now, calcLET is set in pln.bioModels

    pln_HIT.(fn{j}).numOfFractions       = 30;
    pln_HIT.(fn{j}).propStf.gantryAngles = 90; % for the box phantom
    pln_HIT.(fn{j}).propStf.couchAngles   = zeros(numel(pln_HIT.(fn{j}).propStf.gantryAngles),1);

    pln_HIT.(fn{j}).propStf.bixelWidth    = 3; % 3 for MCF and 5 for HIT

    pln_HIT.(fn{j}).propStf.numOfBeams    = numel(pln_HIT.(fn{j}).propStf.gantryAngles);

    pln_HIT.(fn{j}).propStf.isoCenter     = ones(pln_HIT.(fn{j}).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_models.(fn{j}),ct,0);
    
    pln_HIT.(fn{j}).propStf.longitudinalSpotSpacing = 1;
    
    pln_HIT.(fn{j}).propOpt.runDAO        = 0;
    pln_HIT.(fn{j}).propSeq.runSequencing = 0;

    % dose calculation settings
    pln_HIT.(fn{j}).propDoseCalc.doseGrid.resolution.x = 3;
    pln_HIT.(fn{j}).propDoseCalc.doseGrid.resolution.y = 3;
    pln_HIT.(fn{j}).propDoseCalc.doseGrid.resolution.z = 3;

    pln_HIT.(fn{j}).multScen = matRad_multScen(ct, 'nomScen');
end

%% stf generation

% for all bixels

for j = 1: numel(fn)
    stf_HIT.(fn{j}) = matRad_generateStf(ct,cst_models.(fn{j}),pln_HIT.(fn{j}));
end

% for a single bixel
% for j = 1:numel(fn)
%     stf_carbon_HIT.(fn{j}) = matRad_generateSingleBixelStf(ct,cst_models.(fn{j}),pln_HIT.(fn{j}));
% end

%% Pln RBE model

for j = 1:numel(fn)

    if strcmp(fn{j},'LEMI')
        pln_HIT.(fn{j}).bioModel = matRad_bioModel(pln_HIT.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln_HIT.(fn{j}).bioModel.quantityTableName      = 'RBEtable_LEMI_AX217_BX615';

        pln_HIT.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln_HIT.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);

    end

    pln_HIT.(fn{j}).bioModel.calcLET = 1;

    pln_HIT.(fn{j}).propOpt.quantityOpt = 'physicalDose';
    pln_HIT.(fn{j}).propDoseCalc.engine = 'HongPB';
end

%% Dose calculation

dij_HIT.LEMI = matRad_calcDoseInfluence(ct,cst_models.(fn{1}),stf_HIT.(fn{1}),pln_HIT.(fn{1}));

%% Optimization

for j = 1: numel(fn)
    for k = 1:size(cst_models.(fn{j}),1)
        if isempty(cst_models.(fn{j}){k,6})
            cst_models.(fn{j}){k,6} = {};
        else
            cst_models.(fn{j}){k,6}{1}.quantity = 'physicalDose';
        end
    end
end

resultGUI_HIT.LEMI = matRad_fluenceOptimization(dij_HIT.LEMI,cst_models.LEMI,pln_HIT.LEMI);

%% Plot the SOBP

f = figure('WindowState','maximized');

% figure settings
t = tiledlayout(2,3,'TileSpacing','tight','Padding','loose');

set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

set(groot,'defaultAxesFontSize',14)
set(groot,'defaultTextFontSize',14)

% defining the colors
blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
purple = [1 0.1840 0.5560];
brown = [0.6350 0.0780 0.1840];
blueD = [0 0 1];
colors = [blue; orange; green; red; purple; brown; blueD];

% use the central line for the total SOBP and the integral "SOBP" for a
% single bixel plan

lineWidth = 2;
markerSize = 15;
yIndex   = 80;
x = dij_HIT.LEMI.ctGrid.x;

% look at the plot physicalDose in depth before cutting the begin and the
% ending
x = x(41:121);
x = x + 120;

% Plot Physical dose (a)
nexttile;
for j = 1:numel(fn)
    profile_PB_d.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).physicalDose(yIndex,:,yIndex))));
    plot(x,profile_PB_d.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    % profile_MCF_d.(fn{j}) = flip((squeeze(resultGUI_MCF_1456.(fn{j}).physicalDose(yIndex,:,yIndex))));
    % plot(x,profile_MCF_d.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

    % for the integral SOBP
    % integralHITx_d.(fn{j}) = flip(squeeze(sum(resultGUI_HIT.(fn{j}).physicalDose(:,:,:), [2,3])));
    % plot(x,integralHITx_d.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on
    % 
    % integralMCFx_d.(fn{j}) = flip(squeeze(sum(resultGUI_MCF.(fn{j}).physicalDose(:,:,:), [2,3])));
    % plot(x,integralMCFx_d.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on
    
end

grid on;
xlabel('')
xticklabels([])
ylabel('$D_{phys}$ [Gy]');
text(0.95,0.95,'(a)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,2.5])
ax = gca();

% Plot biological dose (b)
nexttile;
for j = 1:numel(fn)
    profile_PB_r.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_PB_r.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    % profile_MCF_r.(fn{j}) = flip((squeeze(resultGUI_MCF_1456.(fn{j}).RBExDose(yIndex,:,yIndex))));
    % plot(x,profile_MCF_r.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

    % for the integral SOBP
    % integralHITx_r.(fn{j}) = flip(squeeze(sum(resultGUI_HIT.(fn{j}).RBExDose(:,:,:), [2,3])));
    % plot(x,integralHITx_r.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on
    % 
    % integralMCFx_r.(fn{j}) = flip(squeeze(sum(resultGUI_MCF.(fn{j}).RBExDose(:,:,:), [2,3])));
    % plot(x,integralMCFx_r.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

end

grid on;
xlabel('')
xticklabels([])
ylabel('$D_{bio,\gamma}$ [Gy(RBE)]');
text(0.95,0.95,'(b)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,6])
ax = gca();

% Plot biological dose carbon (c)
%nexttile;
%for j = [1 2 3]
    % profile_PB_zs.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).zStarDoseCube(yIndex,:,yIndex))));
    % plot(x,profile_PB_zs.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on
    % profile_MCF_zs.(fn{j}) = flip((squeeze(resultGUI_MCF.(fn{j}).zStarDoseCube(yIndex,:,yIndex))));
    % plot(x,profile_MCF_zs.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

    % for the integral SOBP
    % integralHITx_z.(fn{j}) = flip(squeeze(sum(resultGUI_carbon_HIT.(fn{j}).RBExDose(:,:,:), [2,3])));
    % plot(x,integralHITx_z.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on
    % 
    % integralMCFx_z.(fn{j}) = flip(squeeze(sum(resultGUI_carbon_MCF.(fn{j}).RBExDose(:,:,:), [2,3])));
    % plot(x,integralMCFx_z.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

%end
% 
% grid on;
% xlabel('')
% xticklabels([])
% ylabel('$D_{bio,c}$ [Gy(RBE)]');
% text(0.95,0.95,'(c)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
% axis([0,250,0,300])
% ax = gca();

% Plot alpha times dose (d)
nexttile;
for j = 1:numel(fn)
    profile_PB_a.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).alphaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_a.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    % profile_MCF_a.(fn{j}) = flip((squeeze(resultGUI_MCF_1456.(fn{j}).alphaDoseCube(yIndex,:,yIndex))));
    % plot(x,profile_MCF_a.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

    % for the integral SOBP
    % integralHITx_a.(fn{j}) = flip(squeeze(sum(resultGUI_HIT.(fn{j}).alphaDoseCube(:,:,:), [2,3])));
    % plot(x,integralHITx_a.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on
    % 
    % integralMCFx_a.(fn{j}) = flip(squeeze(sum(resultGUI_MCF.(fn{j}).alphaDoseCube(:,:,:), [2,3])));
    % plot(x,integralMCFx_a.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

end

grid on;
xlabel('depth [mm]');
ylabel('$\alpha \times$ dose');
text(0.95,0.95,'(d)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,3.5])
ax = gca();

% Plot sqrt Beta times dose (e)
nexttile;
for j = 1:numel(fn)
    profile_PB_b.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).SqrtBetaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_b.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    % profile_MCF_b.(fn{j}) = flip((squeeze(resultGUI_MCF_1456.(fn{j}).SqrtBetaDoseCube(yIndex,:,yIndex))));
    % plot(x,profile_MCF_b.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

    % for the integral SOBP
    % integralHITx_b.(fn{j}) = flip(squeeze(sum(resultGUI_HIT.(fn{j}).SqrtBetaDoseCube(:,:,:), [2,3])));
    % plot(x,integralHITx_b.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on
    % 
    % integralMCFx_b.(fn{j}) = flip(squeeze(sum(resultGUI_MCF.(fn{j}).SqrtBetaDoseCube(:,:,:), [2,3])));
    % plot(x,integralMCFx_b.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

end

grid on;
xlabel('depth [mm]');
ylabel('$\sqrt{\beta} \times$ dose');
text(0.95,0.95,'(e)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,0.6])
ax = gca();

% Plot LET (f)
nexttile;
for j = 1:numel(fn)
    profile_PB_e.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).LET(yIndex,:,yIndex))));
    plot(x,profile_PB_e.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    % profile_MCF_e.(fn{j}) = flip((squeeze(resultGUI_MCF_1456.(fn{j}).LET(yIndex,:,yIndex))));
    % plot(x,profile_MCF_e.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

    % for the integral SOBP
    % integralHITx_e.(fn{j}) = flip(squeeze(sum(resultGUI_HIT.(fn{j}).LET(:,:,:), [2,3])));
    % plot(x,integralHITx_e.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on
    % 
    % integralMCFx_e.(fn{j}) = flip(squeeze(sum(resultGUI_MCF.(fn{j}).LET(:,:,:), [2,3])));
    % plot(x,integralMCFx_e.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    % hold on

end

grid on;
xlabel('depth [mm]');
ylabel('LET [kev/µm]');
text(0.95,0.95,'(f)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,200])
ax = gca();

% Legend
legend({'LEM-I HIT'}, ...
       'Interpreter','latex', ...
       'Orientation','horizontal')

lgd.Layout.Tile = 'south';

set(gcf,'Position',[100 100 900 500])