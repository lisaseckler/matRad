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
    pln.(fn{j}).machine       = 'MCF_Fluence_updated_withoutRBE_origSingleNewSigma';
    pln.(fn{j}).propDoseCalc.calcLET = 0;

    pln.(fn{j}).numOfFractions       = 30;
    pln.(fn{j}).propStf.gantryAngles = 90; % for the box phantom
    pln.(fn{j}).propStf.couchAngles   = zeros(numel(pln.(fn{j}).propStf.gantryAngles),1);
    pln.(fn{j}).propStf.bixelWidth    = 3;
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
%stf.LEMI = matRad_generateSingleBixelStf(ct,cst_models.LEMI,pln.LEMI);

%% Pln RBE model

for j = 1: numel(fn)

    if strcmp(fn{j},'LEMI')
        pln.(fn{j}).bioModel = matRad_bioModel(pln.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');
        %pln.(fn{j}).bioModel.quantityTableName      = 'RBEtable_LEMI_Scholz06_AX313_BX615';
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

% TabulatedDoseAveragedKernelModel is modified: dE/dx is transposed because
% of MCF_Fluence machine --> denumerator and numerator
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
x = dij_carbon.LEMI.ctGrid.x;
% look at the plot physicalDose in depth before cutting the begin and the
% ending
x = x(41:121);
% x(x < -120) = [];
% x(x > 123) = [];
x = x + 120;

nexttile;
for j = [1 2 3 4]
    profile_PB_d.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).physicalDose(yIndex,:,yIndex))));
    plot(x,profile_PB_d.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_d.(fn{j}) = flip((squeeze(resultGUI_MCF.(fn{j}).physicalDose(yIndex,:,yIndex))));
    plot(x,profile_MCF_d.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
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

for j = [1 2 3 4]
    profile_PB_r.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_PB_r.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_r.(fn{j}) = flip((squeeze(resultGUI_MCF.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_MCF_r.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
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

nexttile;

for j = [1 2 3 4]
    profile_PB_rc.(fn{j}) = flip((squeeze(resultGUI_carbon_HIT.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_PB_rc.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_rc.(fn{j}) = flip((squeeze(resultGUI_carbon_MCF.(fn{j}).RBExDose(yIndex,:,yIndex))));
    plot(x,profile_MCF_rc.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end


grid on;
xlabel('')
xticklabels([])
ylabel('$D_{bio,c}$ [Gy(RBE)]');
text(0.95,0.95,'(c)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,3.5])
ax = gca();

%subplot(2,3,4)
nexttile;
for j = [1 2 3 4]
    profile_PB_a.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).alphaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_a.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_a.(fn{j}) = flip((squeeze(resultGUI_MCF.(fn{j}).alphaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_MCF_a.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('depth [mm]');
ylabel('$\alpha \times$ dose');
text(0.95,0.95,'(d)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,3.5])
ax = gca();

%subplot(2,3,5)
nexttile;
for j = [1 2 3 4]
    profile_PB_b.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).SqrtBetaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_PB_b.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_b.(fn{j}) = flip((squeeze(resultGUI_MCF.(fn{j}).SqrtBetaDoseCube(yIndex,:,yIndex))));
    plot(x,profile_MCF_b.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('depth [mm]');
ylabel('$\sqrt{\beta} \times$ dose');
text(0.95,0.95,'(e)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,0.6])
ax = gca();

%subplot(2,3,6)
nexttile;

for j = [1 2 3 4]
    profile_PB_e.(fn{j}) = flip((squeeze(resultGUI_HIT.(fn{j}).effect(yIndex,:,yIndex))));
    plot(x,profile_PB_e.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    profile_MCF_e.(fn{j}) = flip((squeeze(resultGUI_MCF.(fn{j}).effect(yIndex,:,yIndex))));
    plot(x,profile_MCF_e.(fn{j})(41:121),'--', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end

grid on;
xlabel('depth [mm]');
ylabel('effect');
text(0.95,0.95,'(f)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,240,0,3.5])
ax = gca();

%subplot(2,3,3)

% for j = 1:numel(fn)
%     BED.(fn{j}).BED = resultGUI.(fn{j}).effect ./ cst_models.(fn{j}){k,5}.bioParams.alphaX;
%     profile_PB_bed.(fn{j}) = flip((squeeze(BED.(fn{j}).BED(yIndex,:,yIndex))));
%     plot(x,profile_PB_bed.(fn{j})(41:121),'-', 'DisplayName',num2str(j),'Color', colors(j,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
%     hold on
% end

legend({'LEM-I HIT','LEM-I','mMKM HIT','mMKM','approx. MCF-MKM HIT','approx. MCF-MKM','MCF-MKM HIT','MCF-MKM'}, ...
       'Interpreter','latex', ...
       'Orientation','horizontal')

lgd.Layout.Tile = 'south';

set(gcf,'Position',[100 100 900 500])
%% plot Slice Wrapper

plane = 3;
slice = 80;

LEMI_physD = resultGUI.LEMI.physicalDose;
%LEMIc_physD = resultGUI_RBExDose.LEMI.physicalDose;

% LEMII_physD = resultGUI.LEMII.physicalDose;
% LEMIII_physD = resultGUI.LEMIII.physicalDose;
% MKM_physD = resultGUI.MKM.physicalDose;
mMKM_physD = resultGUI.mMKM.physicalDose;
aMCFMKM_physD = resultGUI.aMCFMKM.physicalDose;
MCFMKM_physD = resultGUI.MCFMKM.physicalDose;

% mMKMc_physD = resultGUI_RBExDose.mMKM.physicalDose;
% aMCFMKMc_physD = resultGUI_RBExDose.aMCFMKM.physicalDose;
% MCFMKMc_physD = resultGUI_RBExDose.MCFMKM.physicalDose;


LEMI_RBExD = resultGUI.LEMI.RBExDose;
% LEMIc_RBExD = resultGUI_RBExDose.LEMI.RBExDose;

% LEMII_RBExD = resultGUI.LEMII.RBExDose;
% LEMIII_RBExD = resultGUI.LEMIII.RBExDose;
% MKM_RBExD = resultGUI.MKM.RBExDose;
mMKM_RBExD = resultGUI.mMKM.RBExDose;
aMCFMKM_RBExD = resultGUI.aMCFMKM.RBExDose;
MCFMKM_RBExD = resultGUI.MCFMKM.RBExDose;

% mMKMc_RBExD = resultGUI_RBExDose.mMKM.RBExDose;
% aMCFMKMc_RBExD = resultGUI_RBExDose.aMCFMKM.RBExDose;
% MCFMKMc_RBExD = resultGUI_RBExDose.MCFMKM.RBExDose;


LEMI_effect = resultGUI.LEMI.effect;
%LEMIc_effect = resultGUI_RBExDose.LEMI.effect;

% LEMII_effect = resultGUI.LEMII.effect;
% LEMIII_effect = resultGUI.LEMIII.effect;
% MKM_effect = resultGUI.MKM.effect;
mMKM_effect = resultGUI.mMKM.effect;
aMCFMKM_effect = resultGUI.aMCFMKM.effect;
MCFMKM_effect = resultGUI.MCFMKM.effect;

% mMKMc_effect = resultGUI_RBExDose.mMKM.effect;
% aMCFMKMc_effect = resultGUI_RBExDose.aMCFMKM.effect;
% MCFMKMc_effect = resultGUI_RBExDose.MCFMKM.effect;

%%
f = figure('WindowState','maximized');

t = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

set(groot,'defaultAxesFontSize',14)
set(groot,'defaultTextFontSize',14)

%subplot(3,4,1)
nexttile;
doseWindow = [0 3];

matRad_plotSliceWrapper(gca,ct,cst,1,LEMI_physD,plane,slice,[],[],colorcube,[],doseWindow);
title('LEM-I')
xlabel('')
xticklabels([])
yticklabels([])
colorbar off
text(0.95,0.95,'(a)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
ylabel('$D_{phys}$')
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
doseWindow = [0 6];

matRad_plotSliceWrapper(gca,ct,cst,1,LEMI_RBExD,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
title('')
colorbar off
text(0.95,0.95,'(e)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
ylabel('$D_{bio,\gamma}$')
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

matRad_plotSliceWrapper(gca,ct,cst,1,LEMI_effect,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
title('')
colorbar off
text(0.95,0.95,'(i)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
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
text(0.95,0.95,'(j)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
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
text(0.95,0.95,'(k)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)

%subplot(3,4,12)
nexttile;
matRad_plotSliceWrapper(gca,ct,cst,1,MCFMKM_effect,plane,slice,[],[],colorcube,[],doseWindow);
xlabel('')
xticklabels([])
yticklabels([])
ylabel('')
title('')
text(0.95,0.95,'(l)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
zoom(2)