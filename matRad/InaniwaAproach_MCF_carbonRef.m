%% Build phantom for Inaniwa2015 reference conditions

ctDim = [90 90 90];

pb = matRad_PhantomBuilder(ctDim,[3 3 3],1);
objective1 = struct(DoseObjectives.matRad_SquaredDeviation(800,60));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));
%pb.addBoxOAR("INSERT",[70 25 25],'offset',[-35 -12 -12],'HU',500)
pb.addBoxTarget("TARGET",[19 19 19],'HU',0,'offset',[16 0 0],'objectives',objective1);
pb.addBoxOAR("OAR",[90,40,40],'HU',0);

[ct, cst] = pb.getctcst();

%% set constraints
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));

%%
% angle 270
modelForComparison = ["LEMI","mMKM","aMCFMKM","MCFMKM"];

for i = 1:size(modelForComparison,2)
    cst_models.(modelForComparison(i)) = cst;
end

% deleting alphaX and betaX in phantom

fn = fieldnames(cst_models);
bioParams = {'alphaX','betaX'};

for j = 1: numel(fn)
    for k = 1:size(cst_models.(fn{j}),1)
        if isfield(cst_models.(fn{j}){k,5},'alphaX') && isfield(cst_models.(fn{j}){k,5},'betaX')
            cst_models.(fn{j}){k,5} = rmfield(cst_models.(fn{j}){k,5},bioParams);
        end
        cst_models.(fn{j}){k,5}.bioParams = struct();
        cst_models.(fn{j}){k,5}.bioParams.cellLine = 'HSG';
        cst_models.(fn{j}){k,5}.bioParams.refRadiation = 'photons'; % carbon if Ztable is used
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
    pln.(fn{j}).propStf.gantryAngles = 270;
    pln.(fn{j}).propStf.couchAngles   = zeros(numel(pln.(fn{j}).propStf.gantryAngles),1);
    pln.(fn{j}).propStf.bixelWidth    = 3;
    pln.(fn{j}).propStf.numOfBeams    = numel(pln.(fn{j}).propStf.gantryAngles);

    pln.(fn{j}).propStf.isoCenter     = ones(pln.(fn{j}).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_models.(fn{j}),ct,0);%
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
for j = 1:numel(fn)
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

dij.MCFMKM = matRad_calcDoseInfluence(ct,cst_models.(fn{4}),stf.(fn{4}),pln.(fn{4}));


%% Optimization

for j = 1: numel(fn)
    for k = 1:size(cst_models.(fn{j}),1)
        if isempty(cst_models.(fn{j}){k,6})
            cst_models.(fn{j}){k,6} = {};
        else
            cst_models.(fn{j}){k,6}{1}.quantity = 'RBExDose';
            %cst_models.(fn{j}){k,6}{2}.quantity = 'physicalDose';
        end
    end
end

resultGUI_RBExD60.(fn{4}) = matRad_fluenceOptimization(dij.(fn{4}),cst_models.(fn{4}),pln.(fn{4}));

%% delete dijs so far
clearvars -regexp ^dij

%% Plot the SOBP

f = figure('WindowState','maximized');
subplot(2,3,1)

blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
purple = [0.4940 0.1840 0.5560];
brown = [0.6350 0.0780 0.1840];
colors = [blue; orange; green; red; purple; brown];

lineWidth = 2;
markerSize = 15;
yIndex   = 45;
x = dij.LEMI.ctGrid.x;
%x = x - 135;

profile_PB_d60.LEMI = ((squeeze(resultGUI_RBExD60.LEMI.physicalDose(yIndex,:,yIndex))));
profile_PB_d60.LEMI(profile_PB_d60.LEMI == 0) = [];

% profile_PB_d60.LEMII = ((squeeze(resultGUI.LEMII.physicalDose(yIndex,:,yIndex))));
% profile_PB_d60.LEMII(profile_PB_d60.LEMII == 0) = [];
% 
% profile_PB_d60.LEMIII = ((squeeze(resultGUI.LEMIII.physicalDose(yIndex,:,yIndex))));
% profile_PB_d60.LEMIII(profile_PB_d60.LEMIII == 0) = [];
% 

profile_PB_d60.aMCFMKM = ((squeeze(resultGUI_RBExD60.aMCFMKM.physicalDose(yIndex,:,yIndex))));
profile_PB_d60.aMCFMKM(profile_PB_d60.aMCFMKM == 0) = [];

profile_PB_d60.MCFMKM = ((squeeze(resultGUI_RBExD60.MCFMKM.physicalDose(yIndex,:,yIndex))));
profile_PB_d60.MCFMKM(profile_PB_d60.MCFMKM == 0) = [];

profile_PB_d60.mMKM = ((squeeze(resultGUI_RBExD60.mMKM.physicalDose(yIndex,:,yIndex))));
profile_PB_d60.mMKM(profile_PB_d60.mMKM == 0) = [];

% profile_PB_d60.MKM = ((squeeze(resultGUI.MKM.physicalDose(yIndex,:,yIndex))));
% profile_PB_d60.MKM(profile_PB_d60.MKM == 0) = [];
% 
% profile_PB_d60.NewmMKM = ((squeeze(resultGUI_RBExD60.NewmMKM.physicalDose(yIndex,:,yIndex))));
% profile_PB_d60.NewmMKM(profile_PB_d60.NewmMKM == 0) = [];


hold on;
plot(x, profile_PB_d60.MCFMKM, '-', 'DisplayName','MCF MKM PB','Color', colors(1,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_d60.mMKM, '-', 'DisplayName','mMKM PB','Color', colors(4,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_d60.aMCFMKM, '-', 'DisplayName','approx. MCF MKM PB','Color', colors(2,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_d60.LEMI, '-', 'DisplayName','LEMI PB','Color', colors(3,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('physicalDose [Gy]','FontSize',20);
axis([0,270,0,1])


ax = gca();


legend('FontSize',20);


subplot(2,3,2)

RBE.LEMI = squeeze(resultGUI_RBExD60.LEMI.RBE(yIndex,:,yIndex));
profile_PB_r60.LEMI = ((squeeze(resultGUI_RBExD60.LEMI.RBExDose(yIndex,:,yIndex))));
profile_newRBExD.LEMI = profile_PB_r60.LEMI ./ RBE.LEMI(61);
profile_newRBExD.LEMI(profile_newRBExD.LEMI == 0) = [];

RBE.aMCFMKM = squeeze(resultGUI_RBExD60.aMCFMKM.RBE(yIndex,:,yIndex));
profile_PB_r60.aMCFMKM = ((squeeze(resultGUI_RBExD60.aMCFMKM.RBExDose(yIndex,:,yIndex))));
profile_newRBExD.aMCFMKM = profile_PB_r60.aMCFMKM ./ RBE.aMCFMKM(61);
profile_newRBExD.aMCFMKM(profile_newRBExD.aMCFMKM == 0) = [];

RBE.MCFMKM = squeeze(resultGUI_RBExD60.MCFMKM.RBE(yIndex,:,yIndex));
profile_PB_r60.MCFMKM = ((squeeze(resultGUI_RBExD60.MCFMKM.RBExDose(yIndex,:,yIndex))));
profile_newRBExD.MCFMKM = profile_PB_r60.MCFMKM ./ RBE.MCFMKM(61);
profile_newRBExD.MCFMKM(profile_newRBExD.MCFMKM == 0) = [];

RBE.mMKM = squeeze(resultGUI_RBExD60.mMKM.RBE(yIndex,:,yIndex));
profile_PB_r60.mMKM = ((squeeze(resultGUI_RBExD60.mMKM.RBExDose(yIndex,:,yIndex))));
profile_newRBExD.mMKM = profile_PB_r60.mMKM ./ RBE.mMKM(61);
profile_newRBExD.mMKM(profile_newRBExD.mMKM == 0) = [];


hold on;
plot(x, profile_newRBExD.MCFMKM, '-', 'DisplayName','MCF MKM PB','Color', colors(1,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newRBExD.mMKM, '-', 'DisplayName','mMKM PB','Color', colors(4,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newRBExD.aMCFMKM, '-', 'DisplayName','approx. MCF MKM PB','Color', colors(2,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newRBExD.LEMI, '-', 'DisplayName','LEMI PB','Color', colors(3,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('biologicalDose [Gy(RBE)]','FontSize',20);
axis([0,270,0,3])

ax = gca();

legend('FontSize',20);

subplot(2,3,3)

profile_PB_e.LEMI = ((squeeze(resultGUI_RBExD60.LEMI.effect(yIndex,:,yIndex))));
profile_PB_e.LEMI(profile_PB_e.LEMI == 0) = [];

profile_PB_e.aMCFMKM = ((squeeze(resultGUI_RBExD60.aMCFMKM.effect(yIndex,:,yIndex))));
profile_PB_e.aMCFMKM(profile_PB_e.aMCFMKM == 0) = [];

profile_PB_e.MCFMKM = ((squeeze(resultGUI_RBExD60.MCFMKM.effect(yIndex,:,yIndex))));
profile_PB_e.MCFMKM(profile_PB_e.MCFMKM == 0) = [];

profile_PB_e.mMKM = ((squeeze(resultGUI_RBExD60.mMKM.effect(yIndex,:,yIndex))));
profile_PB_e.mMKM(profile_PB_e.mMKM == 0) = [];

hold on;
plot(x, profile_PB_e.MCFMKM, '-', 'DisplayName','MCF MKM PB','Color', colors(1,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_e.mMKM, '-', 'DisplayName','mMKM PB','Color', colors(4,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_e.aMCFMKM, '-', 'DisplayName','approx. MCF MKM PB','Color', colors(2,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_e.LEMI, '-', 'DisplayName','LEMI PB','Color', colors(3,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('effect','FontSize',20);
axis([0,270,0,1.5])
ax = gca();

legend('FontSize',20);

subplot(2,3,4)

profile_PB_a.LEMI = ((squeeze(resultGUI_RBExD60.LEMI.alphaDoseCube(yIndex,:,yIndex))));
profile_newAlpha.LEMI = profile_PB_a.LEMI ./ RBE.LEMI(61);
profile_newAlpha.LEMI(profile_newAlpha.LEMI == 0) = [];

profile_PB_a.aMCFMKM = ((squeeze(resultGUI_RBExD60.aMCFMKM.alphaDoseCube(yIndex,:,yIndex))));
profile_newAlpha.aMCFMKM = profile_PB_a.aMCFMKM ./ RBE.aMCFMKM(61);
profile_newAlpha.aMCFMKM(profile_newAlpha.aMCFMKM == 0) = [];

profile_PB_a.MCFMKM = ((squeeze(resultGUI_RBExD60.MCFMKM.alphaDoseCube(yIndex,:,yIndex))));
profile_newAlpha.MCFMKM = profile_PB_a.MCFMKM ./ RBE.MCFMKM(61);
profile_newAlpha.MCFMKM(profile_newAlpha.MCFMKM == 0) = [];

profile_PB_a.mMKM = ((squeeze(resultGUI_RBExD60.mMKM.alphaDoseCube(yIndex,:,yIndex))));
profile_newAlpha.mMKM = profile_PB_a.mMKM ./ RBE.mMKM(61);
profile_newAlpha.mMKM(profile_newAlpha.mMKM == 0) = [];


hold on;
plot(x, profile_newAlpha.MCFMKM, '-', 'DisplayName','MCF MKM PB','Color', colors(1,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newAlpha.mMKM, '-', 'DisplayName','mMKM PB','Color', colors(4,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newAlpha.aMCFMKM, '-', 'DisplayName','approx. MCF MKM PB','Color', colors(2,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newAlpha.LEMI, '-', 'DisplayName','LEMI PB','Color', colors(3,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('alphaDose','FontSize',20);
axis([0,270,0,1.5])

ax = gca();

legend('FontSize',20);

subplot(2,3,5)

profile_PB_b.MCFMKM = ((squeeze(resultGUI_RBExD60.MCFMKM.SqrtBetaDoseCube(yIndex,:,yIndex))));
profile_newSqrtBeta.MCFMKM = profile_PB_b.MCFMKM ./ RBE.MCFMKM(61);
profile_newSqrtBeta.MCFMKM(profile_newSqrtBeta.MCFMKM == 0) = [];

profile_PB_b.mMKM = ((squeeze(resultGUI_RBExD60.mMKM.SqrtBetaDoseCube(yIndex,:,yIndex))));
profile_newSqrtBeta.mMKM = profile_PB_b.mMKM ./ RBE.mMKM(61);
profile_newSqrtBeta.mMKM(profile_newSqrtBeta.mMKM == 0) = [];

profile_PB_b.aMCFMKM = ((squeeze(resultGUI_RBExD60.aMCFMKM.SqrtBetaDoseCube(yIndex,:,yIndex))));
profile_newSqrtBeta.aMCFMKM = profile_PB_b.aMCFMKM ./ RBE.aMCFMKM(61);
profile_newSqrtBeta.aMCFMKM(profile_newSqrtBeta.aMCFMKM == 0) = [];

profile_PB_b.LEMI = ((squeeze(resultGUI_RBExD60.LEMI.SqrtBetaDoseCube(yIndex,:,yIndex))));
profile_newSqrtBeta.LEMI = profile_PB_b.LEMI ./ RBE.LEMI(61);
profile_newSqrtBeta.LEMI(profile_newSqrtBeta.LEMI == 0) = [];

hold on;
plot(x, profile_newSqrtBeta.MCFMKM, '-', 'DisplayName','MCF MKM PB','Color', colors(1,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newSqrtBeta.mMKM, '-', 'DisplayName','mMKM PB','Color', colors(4,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newSqrtBeta.aMCFMKM, '-', 'DisplayName','approx. MCF MKM PB','Color', colors(2,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_newSqrtBeta.LEMI, '-', 'DisplayName','LEMI PB','Color', colors(3,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('SqrtBetaDose','FontSize',20);
axis([0,270,0,0.3])
ax = gca();
legend('FontSize',20);
%%
subplot(2,3,6)
% for j = 1:4
%     BED.(fn{j}).BED = resultGUI_RBExD60.(fn{j}).effect ./ cst_models.(fn{j}){k,5}.bioParams.alphaX;
% end

BED.NewmMKM.BED = resultGUI_RBExD60.NewmMKM.effect ./ cst_models.mMKM{k,5}.bioParams.alphaX;

profile_PB_R.MCFMKM = ((squeeze(BED.MCFMKM.BED(yIndex,:,yIndex))));
% profile_PB_R.MCFMKM = profile_PB_R.MCFMKM(1:35);

profile_PB_R.mMKM = ((squeeze(BED.mMKM.BED(yIndex,:,yIndex))));
%profile_PB_R.mMKM = profile_PB_R.mMKM(1:35);

profile_PB_R.NewmMKM = ((squeeze(BED.NewmMKM.BED(yIndex,:,yIndex))));

profile_PB_R.aMCFMKM = ((squeeze(BED.aMCFMKM.BED(yIndex,:,yIndex))));
%profile_PB_R.aMCFMKM = profile_PB_R.aMCFMKM(1:35);

profile_PB_R.LEMI = ((squeeze(BED.LEMI.BED(yIndex,:,yIndex))));
%profile_PB_R.LEMI = profile_PB_R.LEMI(1:35);

hold on;
plot(x, profile_PB_R.MCFMKM, '-', 'DisplayName','MCF MKM PB','Color', colors(1,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_R.mMKM, '-', 'DisplayName','mMKM PB','Color', colors(4,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_R.NewmMKM, '-', 'DisplayName','new mMKM PB','Color', colors(5,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_R.aMCFMKM, '-', 'DisplayName','approx. MCF MKM PB','Color', colors(2,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);
hold on;
plot(x, profile_PB_R.LEMI, '-', 'DisplayName','LEMI PB','Color', colors(3,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize);

grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('BED [Gy]','FontSize',20);
axis([0,270,0,5])

ax = gca();

legend('FontSize',20);

%% alpha and beta Ref for carbon reference

for j = 1:4
    alphaRef.(fn{j}) = profile_PB_a.(fn{j})(61) / profile_PB_d60.(fn{j})(61);
    betaRef.(fn{j}) = (profile_PB_b.(fn{j})(61) / profile_PB_d60.(fn{j})(61))^2;
end

%%
alphaRef.NewmMKM = profile_PB_a.NewmMKM(61) / profile_PB_d60.NewmMKM(61);
betaRef.NewmMKM = profile_PB_b.NewmMKM(61) / profile_PB_d60.NewmMKM(61);

%% saving
save("BioParamMCF_CarbonRef_18032026.mat","alphaRef","betaRef")