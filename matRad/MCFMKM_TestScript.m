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
modelForComparison = ["MCFMKM"];

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
dij.MCFMKM = matRad_calcDoseInfluence(ct,cst_models.(fn{1}),stf.(fn{1}),pln.(fn{1}));

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

resultGUI_RBExD60.MCFMKM = matRad_fluenceOptimization(dij.MCFMKM,cst_models.MCFMKM,pln.MCFMKM);