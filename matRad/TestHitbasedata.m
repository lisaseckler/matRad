%% test HIT base data

load('BOXPHANTOM.mat');

modelForComparison = ["LEMI"];

for i = 1:size(modelForComparison,2)
    cst_same.(modelForComparison(i)) = cst;
end

% deleting alphaX and betaX in phantom

fn = fieldnames(cst_same);
bioParams = {'alphaX','betaX'};

for j = 1: numel(fn)
    for k = 1:size(cst_same.(fn{j}),1)
        if isfield(cst_same.(fn{j}){k,5},'alphaX') && isfield(cst_same.(fn{j}){k,5},'betaX')
            cst_same.(fn{j}){k,5} = rmfield(cst_same.(fn{j}){k,5},bioParams);
        end
        cst_same.(fn{j}){k,5}.bioParams = struct();
        cst_same.(fn{j}){k,5}.bioParams.cellLine = 'HSG';
        cst_same.(fn{j}){k,5}.bioParams.refRadiation = 'photons'; % carbon if Ztable is used
        cst_same.(fn{j}){k,5}.bioParams.alphaX = 0.313;
        cst_same.(fn{j}){k,5}.bioParams.betaX = 0.0615;
    end
end


%% Pln Definition

for j = 1: numel(fn)

    pln_same.(fn{j}).radiationMode = 'carbon';
    pln_same.(fn{j}).machine       = 'HITfluenceSpectra_updatedBig';
    %pln_prostate.(fn{j}).machine       = 'Generic_clusterDose_prestep_update';
    pln_same.(fn{j}).propDoseCalc.calcLET = 0;

    pln_same.(fn{j}).numOfFractions       = 30;
    pln_same.(fn{j}).propStf.gantryAngles = 90;%[90 270];
    pln_same.(fn{j}).propStf.couchAngles   = zeros(numel(pln_same.(fn{j}).propStf.gantryAngles),1);
    pln_same.(fn{j}).propStf.bixelWidth    = 5;
    pln_same.(fn{j}).propStf.numOfBeams    = numel(pln_same.(fn{j}).propStf.gantryAngles);

    pln_same.(fn{j}).propStf.isoCenter     = ones(pln_same.(fn{j}).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_same.(fn{j}),ct,0);
    pln_same.(fn{j}).propOpt.runDAO        = 0;
    pln_same.(fn{j}).propSeq.runSequencing = 0;

    % dose calculation settings
    pln_same.(fn{j}).propDoseCalc.doseGrid.resolution.x = 3;
    pln_same.(fn{j}).propDoseCalc.doseGrid.resolution.y = 3;
    pln_same.(fn{j}).propDoseCalc.doseGrid.resolution.z = 3;

    pln_same.(fn{j}).multScen = matRad_multScen(ct, 'nomScen');
end


%% stf generation
for j = 1: numel(fn)
    stf_same.(fn{j}) = matRad_generateStf(ct,cst_same.(fn{j}),pln_same.(fn{j}));
end

for j = 1: numel(fn)
    stf_oldMachine.(fn{j}) = matRad_generateStf(ct,cst_same.(fn{j}),pln_oldMachine.(fn{j}));
end
%stf = matRad_generateSingleBixelStf(ct,cst,pln);

%% Pln RBE model

for j = 1: numel(fn)
    pln_same.(fn{j}).bioModel = matRad_bioModel(pln_same.(fn{j}).radiationMode, 'doseAveragedTabulatedAlphaBeta');

    if strcmp(fn{j},'LEMI')
       pln_same.(fn{j}).bioModel.quantityTableName      = 'RBEtable_LEMI_Scholz06_AX313_BX615';
    end
    
    pln_same.(fn{j}).bioModel.stoppingPowerTableName = 'SPtable';

    selectZ = [1,2,3,4,5,6];

    selectA = [1,4,7,9,11,12];

    pln_same.(fn{j}).bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);

    pln_same.(fn{j}).propOpt.quantityOpt = 'physicalDose';
    pln_same.(fn{j}).propDoseCalc.engine = 'HongPB';
end
%%
pln_oldMachine.LEMI.bioModel = matRad_bioModel(pln_oldMachine.LEMI.radiationMode,'LEM');

%% Dose calculation

% everything is based on the LEMI optimization (weight vector)
dij_old = matRad_calcDoseInfluence(ct,cst_same.(fn{1}),stf_oldMachine.(fn{1}),pln_oldMachine.(fn{1}));

%% Optimization

for j = 1: numel(fn)
    for k = 1:size(cst_same.(fn{j}),1)
        if isempty(cst_same.(fn{j}){k,6})
            cst_same.(fn{j}){k,6} = {};
        else
            cst_same.(fn{j}){k,6}{1}.quantity = 'physicalDose';
        end
    end
end

resultGUI_oldMachine.LEMI2 = matRad_fluenceOptimization(dij_old,cst_same.LEMI,pln_oldMachine.LEMI);

%% Recalculation
pln_oldMachine.LEMI = pln_same.LEMI;
pln_oldMachine.LEMI.machine = 'HITfixedBL';

resultGUI_oldMachine.LEMI = matRad_calcDoseForward(ct,cst_same.LEMI,stf_same.LEMI,pln_oldMachine.LEMI,resultGUI_HITspectra.LEMI.w);

%% Gamma Index

resultGUI_append2 = matRad_appendResultGUI(resultGUI_HITspectra.LEMI,resultGUI_oldMachine.LEMI2,true,pln_same.LEMI.propDoseCalc.engine);
compDose2 = matRad_compareDose(resultGUI_append2.physicalDose, resultGUI_append2.(['physicalDose_' pln_same.LEMI.propDoseCalc.engine]), ct, cst_same.LEMI, [1, 1, 0] , 'off', pln_same.LEMI, [3, 3], 3, 'global');

%% SOBPs

f = figure('WindowState','maximized');
subplot(1,2,1)

blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
purple = [0.4940 0.1840 0.5560];
brown = [0.6350 0.0780 0.1840];
colors = [blue; orange; green; red; purple; brown];

lineWidth = 2;
markerSize = 15;
yIndex   = 80;
x = dij_same.ctGrid.x;
% x(x < -115) = [];
% x(x > 3) = [];
% x = x + 115;

profile_PB_d_LEMI = flip((squeeze(resultGUI_HITspectra.LEMI.physicalDose(yIndex,:,yIndex))));
profile_PB_d_LEMI2 = flip((squeeze(resultGUI_oldMachine.LEMI2.physicalDose(yIndex,:,yIndex))));

hold on;
plot(x, profile_PB_d_LEMI, '-', 'DisplayName','LEMI spectra','Color', colors(3,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on
plot(x, profile_PB_d_LEMI2, '--', 'DisplayName','LEMI old','Color', colors(2,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('physicalDose [Gy]','FontSize',20);
axis([-240,240,0,2.5])


ax = gca();


legend('FontSize',20);


subplot(1,2,2)

profile_PB_r_LEMI = flip((squeeze(resultGUI_HITspectra.LEMI.RBExDose(yIndex,:,yIndex))));
profile_PB_r_LEMI2 = flip((squeeze(resultGUI_oldMachine.LEMI2.RBExDose(yIndex,:,yIndex))));

hold on;
plot(x, profile_PB_r_LEMI, '-', 'DisplayName','LEMI spectra','Color', colors(3,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_r_LEMI2, '--', 'DisplayName','LEMI old','Color', colors(2,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('biologicalDose [Gy(RBE)]','FontSize',20);
axis([-240,240,0,6])

ax = gca();

legend('FontSize',20);
