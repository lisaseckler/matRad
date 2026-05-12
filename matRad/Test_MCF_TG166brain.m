%% Master Thesis: TG166 brain comparisons based on --> photon reference
%matRad_rc;
clear
matRad_cfg = MatRad_Config.instance();

load("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\userdata\patients\TG166_brain_ct3mm.mat")

%% set objectives
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,9));
%cst{1,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));

cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(1000,18));
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,16));

cst{6,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(200,9));
cst{14,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,9));

%% define models to compare

modelForComparison = ["LEMI","mMKM","aMCFMKM","MCFMKM"];

for i = 1:size(modelForComparison,2)
    cst_models_MCF_carbon.(modelForComparison(i)) = cst;
end

% deleting alphaX and betaX in phantom and setting model specific reference
% parameters

fn = fieldnames(cst_models_MCF_carbon);
bioParams = {'alphaX','betaX'};

for j = 1: numel(fn)
    for k = 1:size(cst_models_MCF_carbon.(fn{j}),1)
        if isfield(cst_models_MCF_carbon.(fn{j}){k,5},'alphaX') && isfield(cst_models_MCF_carbon.(fn{j}){k,5},'betaX')
            cst_models_MCF_carbon.(fn{j}){k,5} = rmfield(cst_models_MCF_carbon.(fn{j}){k,5},bioParams);
        end
        cst_models_MCF_carbon.(fn{j}){k,5}.bioParams = struct();
        cst_models_MCF_carbon.(fn{j}){k,5}.bioParams.cellLine = 'HSG';
        cst_models_MCF_carbon.(fn{j}){k,5}.bioParams.refRadiation = 'photons'; 
        cst_models_MCF_carbon.(fn{j}){k,5}.bioParams.alphaX = 0.217;
        cst_models_MCF_carbon.(fn{j}){k,5}.bioParams.betaX = 0.0615;

    end
end

%% Pln Definition

for j = 1: numel(fn)

    pln.(fn{j}).radiationMode = 'carbon';
    pln.(fn{j}).machine       = 'MCF_Fluence_updated_withoutRBE_origSingleNewSigma';

    pln.(fn{j}).propDoseCalc.calcLET = 0;

    pln.(fn{j}).numOfFractions       = 3;
    pln.(fn{j}).propStf.gantryAngles = [-60 60]; %[-85 -45 45 85];
    pln.(fn{j}).propStf.couchAngles   = zeros(numel(pln.(fn{j}).propStf.gantryAngles),1);
    pln.(fn{j}).propStf.bixelWidth    = 3;
    pln.(fn{j}).propStf.numOfBeams    = numel(pln.(fn{j}).propStf.gantryAngles);

    pln.(fn{j}).propStf.isoCenter     = ones(pln.(fn{j}).propStf.numOfBeams,1) * matRad_getIsoCenter(cst_models_MCF_carbon.(fn{j}),ct,0);
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
    stf_photon.(fn{j}) = matRad_generateStf(ct,cst_models_MCF_carbon.(fn{j}),pln.(fn{j}));
end

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

        %pln.(fn{j}).bioModel.quantityTableName = 'RBEtable_mMKM_Inaniwa10_A0B0_AX3312_BX593';
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

dij_photon.LEMI = matRad_calcDoseInfluence(ct,cst_models_MCF_carbon.(fn{1}),stf_photon.(fn{1}),pln.(fn{1}));

%% Optimization

% set optimization quantity
% for the box phantom the quantity should be physical dose, so the dose
% stays the same and the RBE can be analyzed

for j = 1: numel(fn)
    for k = 1:size(cst_models_MCF_carbon.(fn{j}),1)
        if isempty(cst_models_MCF_carbon.(fn{j}){k,6})
            cst_models_MCF_carbon.(fn{j}){k,6} = {};
        else
            cst_models_MCF_carbon.(fn{j}){k,6}{1}.quantity = 'RBExDose';
        end
    end
end

resultGUI_photon.LEMI = matRad_fluenceOptimization(dij_photon.LEMI,cst_models_MCF_carbon.LEMI,pln.LEMI);


%% Recalculation
% same weights for the same dose for every biological model

for j = 2:numel(fn)
    resultGUI_photon.(fn{j}) = matRad_calcDoseForward(ct,cst_models_MCF_carbon.(fn{j}),stf_photon.(fn{j}),pln.(fn{j}),resultGUI_photon.LEMI.w);
end
