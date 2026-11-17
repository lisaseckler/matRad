%matRad_rc;

% Compare Ztables
clear
matRad_cfg = MatRad_Config.instance();

load("BOXPHANTOM.mat")

%% set objectives
% use standard objectives
% penalty from tumor set to 1000

%% define models to compare
cst{2,6}{1,1}.penalty = 1500;

% deleting alphaX and betaX in phantom and setting model specific reference
% parameters

bioParams = {'alphaX','betaX'};

for k = 1:size(cst,1)
    if isfield(cst{k,5},'alphaX') && isfield(cst{k,5},'betaX')
        cst{k,5} = rmfield(cst{k,5},bioParams);
    end
    cst{k,5}.bioParams = struct();
    cst{k,5}.bioParams.cellLine = 'HSG';
    cst{k,5}.bioParams.refRadiation = 'photons';
    cst{k,5}.bioParams.alphaX = 0.172;
    cst{k,5}.bioParams.alphaR = 0.3312;
    cst{k,5}.bioParams.betaX = 0.0615;
    cst{k,5}.bioParams.betaR = 0.0615;
end

%% Pln Definition

pln.radiationMode = 'carbon';
pln.machine       = 'MCF_doubleGauss_ridgeFilter_sigmaAveraged';
%pln.machine        = 'HITfluenceSpectra_withEnergySpectrum';
pln.propDoseCalc.calcLET = 0;

pln.numOfFractions       = 30;
pln.propStf.gantryAngles = 270; %90; % for the box phantom
pln.propStf.couchAngles   = zeros(numel(pln.propStf.gantryAngles),1);
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);

pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propStf.longitudinalSpotSpacing = 1;
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3;
pln.propDoseCalc.doseGrid.resolution.y = 3;
pln.propDoseCalc.doseGrid.resolution.z = 3;

pln.multScen = matRad_multScen(ct, 'nomScen');

%% stf generation

stf = matRad_generateStf(ct,cst,pln);

%% z* table
% if a z* table is used the table must be called like this
pln.bioModel = matRad_bioModel(pln.radiationMode, 'doseAveragedTabulatedAlphaBeta');
pln.bioModel.ZstarTableName = 'Ztable_mMKM_Shannon';%'Ztable_mMKM_zs_Inaniwa10_A0172_B0615';
pln.bioModel.stoppingPowerTableName = 'SPtable';

selectZ = [1,2,3,4,5,6];

selectA = [1,4,7,9,11,12];

pln.bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);
pln.propOpt.quantityOpt = 'physicalDose';
pln.bioModel.calcLET = 1;
pln.propDoseCalc.engine = 'HongPB';

%% Dose calculation

dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Optimization

% set optimization quantity
% for the box phantom the quantity should be physical dose, so the dose
% stays the same and the RBE can be analyzed


for k = 1:size(cst,1)
    if isempty(cst{k,6})
        cst{k,6} = {};
    else
        cst{k,6}{1}.quantity = 'physicalDose';
    end
end

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the SOBP

f = figure('WindowState','maximized');

t = tiledlayout(2,2,'TileSpacing','tight','Padding','loose');

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
x = dij.ctGrid.x;
% look at the plot physicalDose in depth before cutting the begin and the
% ending
x = x(41:121);
% x(x < -120) = [];
% x(x > 123) = [];
x = x + 120;

nexttile;

profile_me_d = flip((squeeze(resultGUI_me.physicalDose(yIndex,:,yIndex))));
plot(x,profile_me_d(41:121),'-', 'DisplayName','me','Color', [0 0 1], 'LineWidth', lineWidth, 'MarkerSize',markerSize)
hold on
profile_Shannon_d = flip((squeeze(resultGUI.physicalDose(yIndex,:,yIndex))));
plot(x,profile_Shannon_d(41:121),'--', 'DisplayName','Shannon','Color', [0.9,0.6,0], 'LineWidth', lineWidth, 'MarkerSize',markerSize)
hold on

grid on;
% xlabel('depth [mm]', 'FontSize',20);
xlabel('')
xticklabels([])
ylabel('$D_{phys}$ [Gy]');
text(0.95,0.95,'(a)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,2.5])
ax = gca();

nexttile;

profile_me_r = flip((squeeze(resultGUI_me.RBExDose(yIndex,:,yIndex))));
plot(x,profile_me_r(41:121),'-', 'DisplayName','me','Color', [0 0 1], 'LineWidth', lineWidth, 'MarkerSize',markerSize)
hold on
profile_Shannon_r = flip((squeeze(resultGUI.RBExDose(yIndex,:,yIndex))));
plot(x,profile_Shannon_r(41:121),'--', 'DisplayName','Shannon','Color', [0.9,0.6,0], 'LineWidth', lineWidth, 'MarkerSize',markerSize)
hold on

grid on;
%xlabel('depth [mm]', 'FontSize',20);
xlabel('')
xticklabels([])
ylabel('$D_{bio,\gamma}$ [Gy(RBE)]');
text(0.95,0.95,'(b)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,6])
ax = gca();

nexttile;

profile_me_a = flip((squeeze(resultGUI_me.alphaDoseCube(yIndex,:,yIndex))));
plot(x,profile_me_a(41:121),'-', 'DisplayName','me','Color', [0 0 1], 'LineWidth', lineWidth, 'MarkerSize',markerSize)
hold on
profile_Shannon_a = flip((squeeze(resultGUI.alphaDoseCube(yIndex,:,yIndex))));
plot(x,profile_Shannon_a(41:121),'--', 'DisplayName','Shannon','Color', [0.9,0.6,0], 'LineWidth', lineWidth, 'MarkerSize',markerSize)
hold on

grid on;
xlabel('depth [mm]');
ylabel('$\alpha \times$ dose');
text(0.95,0.95,'(d)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,3.5])
ax = gca();

nexttile;

profile_me_b = flip((squeeze(resultGUI_me.SqrtBetaDoseCube(yIndex,:,yIndex))));
plot(x,profile_me_b(41:121),'-', 'DisplayName','me','Color', [0 0 1], 'LineWidth', lineWidth, 'MarkerSize',markerSize)
hold on
profile_Shannon_b = flip((squeeze(resultGUI.SqrtBetaDoseCube(yIndex,:,yIndex))));
plot(x,profile_Shannon_b(41:121),'--', 'DisplayName','Shannon','Color', [0.9,0.6,0], 'LineWidth', lineWidth, 'MarkerSize',markerSize)
hold on

grid on;
xlabel('depth [mm]');
ylabel('$\sqrt{\beta} \times$ dose');
text(0.95,0.95,'(e)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
axis([0,250,0,0.6])
ax = gca();

% Legend
hold on
% --- Dummy handles for MODEL (color) ---
me = plot(nan, nan, '-', 'Color', [0 0 1], 'LineWidth', lineWidth); % j=1, blue
Shannon = plot(nan, nan, '-', 'Color', [0.9,0.6,0], 'LineWidth', lineWidth); % j=2, orange

lgd = legend(...
    {'me','Shannon'}, ...
    'Interpreter','latex', ...
    'Orientation','horizontal');

lgd.Layout.Tile = 'south';

set(gcf,'Position',[100 100 900 500])