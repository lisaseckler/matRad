%% MCF: box phantom comparisons based on --> photon reference
%matRad_rc;
% clear
matRad_cfg = MatRad_Config.instance();

load("BOXPHANTOM.mat")

%% set objectives
% use standard objectives
% penalty from tumor set to 1000

%% define models to compare

modelForComparison = ["LEMI","mMKM","aMCFMKM","MCFMKM"];
cst{2,6}{1,1}.penalty = 1500;

for i = 1:size(modelForComparison,2)
    cst_models.(modelForComparison(i)) = cst;
end

% deleting alphaX and betaX in phantom and setting model specific reference
% parameters

fn = fieldnames(cst_models);
bioParams = {'alphaX','betaX'};

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

for j = 1: numel(fn)

    pln.(fn{j}).radiationMode = 'carbon';
    pln.(fn{j}).machine       = 'MCF_Fluence_updated_withoutRBE_singleNewSigma';
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

for j = 1: numel(fn)

    if strcmp(fn{j},'LEMI')
        pln.(fn{j}).bioModel = matRad_bioModel(pln.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln.(fn{j}).bioModel.quantityTableName      = 'RBEtable_LEMI_AX217_BX615';
        
        pln.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);


    end

    if strcmp(fn{j},'mMKM')
        pln.(fn{j}).bioModel = matRad_bioModel(pln.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln.(fn{j}).bioModel.quantityTableName = 'RBEtable_mMKM_photon_A0_AX217_BX615';

        pln.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

        selectZ = [1,2,3,4,5,6];

        selectA = [1,4,7,9,11,12];

        pln.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);


    end

    if strcmp(fn{j},'aMCFMKM')
        pln.(fn{j}).bioModel = matRad_bioModel(pln.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        pln.(fn{j}).bioModel.quantityTableName = 'RBEtable_MCFMKM_photon_AX0217_BX00615';
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



    pln.(fn{j}).propOpt.quantityOpt = 'physicalDose';
    pln.(fn{j}).propDoseCalc.engine = 'HongPB';
end

%% Dose calculation

% everything is based on the LEMI optimization (weight vector)
% only calculate the dij once for LEMI

dij.LEMI = matRad_calcDoseInfluence(ct,cst_models.(fn{1}),stf.(fn{1}),pln.(fn{1}));

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

resultGUI.LEMI = matRad_fluenceOptimization(dij.LEMI,cst_models.LEMI,pln.LEMI);

%% Recalculation
% same weights for the same dose for every biological model

for j = 2:numel(fn)
    resultGUI.(fn{j}) = matRad_calcDoseForward(ct,cst_models.(fn{j}),stf.(fn{j}),pln.(fn{j}),resultGUI.LEMI.w);
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
x = dij.mMKM.ctGrid.x;
% look at the plot physicalDose in depth before cutting the begin and the
% ending
x = x(41:121);
% x(x < -120) = [];
% x(x > 123) = [];
x = x + 120;

nexttile;
for j = [1 2]
    profile_PB_d.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).physicalDose(yIndex,:,yIndex))));
    plot(x,profile_PB_d.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
% xlabel('depth [mm]', 'FontSize',20);
xlabel('')
xticklabels([])
ylabel('$D_{phys}$ [Gy]');
text(0.95,0.95,'(a)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,2.5])
ax = gca();

%subplot(2,3,2)
nexttile;

for j = [1 3]
    profile_PB_r.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_PB_r.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
%xlabel('depth [mm]', 'FontSize',20);
xlabel('')
xticklabels([])
ylabel('$D_{bio,\gamma}$ [Gy(RBE)]');
text(0.95,0.95,'(b)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,6])
ax = gca();

% nexttile;
% 
% for j = [1 2 3 4]
%     profile_PB_r.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).RBExDose(yIndex,:,yIndex))));
%     plot(x,profile_PB_r.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
%     hold on
% end
% 
% 
% grid on;
% xlabel('')
% xticklabels([])
% ylabel('$D_{bio,c}$ [Gy(RBE)]');
% text(0.95,0.95,'(c)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
% axis([0,240,0,3.5])
% ax = gca();

%subplot(2,3,4)
nexttile;
for j = [1 3]
    profile_PB_a.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).alphaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_a.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('')
xticklabels([])
ylabel('$\alpha \times$ dose');
text(0.95,0.95,'(b)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,3.5])
ax = gca();

%subplot(2,3,5)
nexttile;
for j = [1 3]
    profile_PB_b.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).SqrtBetaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_b.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('depth [mm]');
ylabel('$\sqrt{\beta} \times$ dose');
text(0.95,0.95,'(c)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,0.6])
ax = gca();

%subplot(2,3,6)
nexttile;

for j = [1 3]
    profile_PB_e.(fn{j}) = flip((squeeze(resultGUI.(fn{j}).effect(yIndex,:,yIndex))));
    plot(x,profile_PB_e.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('depth [mm]');
ylabel('effect');
text(0.95,0.95,'(d)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,4])
ax = gca();

%subplot(2,3,3)

% for j = 1:numel(fn)
%     BED.(fn{j}).BED = resultGUI.(fn{j}).effect ./ cst_models.(fn{j}){k,5}.bioParams.alphaX;
%     profile_PB_bed.(fn{j}) = flip((squeeze(BED.(fn{j}).BED(yIndex,:,yIndex))));
%     plot(x,profile_PB_bed.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
%     hold on
% end

legend({'mMKM with alpha-beta table','mMKM with Ztable'}, ...
       'Interpreter','latex', ...
       'Orientation','horizontal')

lgd.Layout.Tile = 'south';

set(gcf,'Position',[100 100 900 500])