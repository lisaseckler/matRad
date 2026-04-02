%% Example: Tabulated BioModels for Carbon Ion Treatment Plan

%% In this example we will show 
% (i) how to set the cst with all the bio parameters that are needed
% (ii) how to set up a carbon ion dose calculation plan with tabulated RBE models
% (iii) how to set up a carbon ion treatment plan with a z*table

%% set matRad runtime configuration
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

%% Patient Data Import
load('BOXPHANTOM.mat');

%% Setting the cst right
% Depending on the cell line that is used for generating an alpha-beta
% table, the cell line needs to be defined in the cst. Additionally, the
% alpha and beta reference values need to be set. To make sure the correct
% reference radiation is used, it's nice to also set the reference
% radiation.
%%
% the for loop is only necessary if the same alpha/beta ratios are used for
% every structure in the cst
for k = 1:size(cst,1)
    cst{k,5}.bioParams = struct(); % create a struct

    % set the cell line: This is important because it will check in the dij calculation
    % if it matches with the cell line that is defined in the alpha-beta table
    cst{k,5}.bioParams.cellLine = 'HSG';

    % the reference radiation does not need to be set, it is just for clarity
    cst{k,5}.bioParams.refRadiation = 'photons';

    % set the alpha and beta values for the chosen reference radiation
    cst{k,5}.bioParams.alphaX = 0.313;
    cst{k,5}.bioParams.betaX = 0.0615;

    % if you choose to use a z* table, then you need alpha0 as well, which is
    % why alphaX and betaX are used as input parameters basically alpha0 and
    % alphaR and betaR describe the reference parameters

    % cst{k,5}.bioParams.alphaX = 0.172;
    % cst{k,5}.bioParams.betaX = 0.0615;
    % cst{k,5}.bioParams.alphaR = 0.3312;
    % cst{k,5}.bioParams.betaR = 0.0615; % betaR does not have to be betaX, it depends on the reference radiation
end

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.
%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are protons, helium or carbon. In this
% example we would like to use carbon ions for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. For tabulated bioModels the fluence need to be scored in the
% base data to average the alpha and beta with the fluence (comparability reasons with TOPAS)

pln.radiationMode   = 'carbon';            
pln.machine         = 'Generic_clusterDose_prestep_update';
pln.numOfFractions  = 30;

pln.propDoseCalc.calcLET = 0;

pln.propStf.gantryAngles = 90;
pln.propStf.couchAngles  = zeros(numel(pln.propStf.gantryAngles),1);
pln.propStf.bixelWidth   = 5;
pln.propStf.numOfBeams   = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propStf.longitudinalSpotSpacing = 1;

pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3;
pln.propDoseCalc.doseGrid.resolution.y = 3;
pln.propDoseCalc.doseGrid.resolution.z = 3;

pln.multScen = matRad_multScen(ct, 'nomScen');

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Set pln bioModel parameters
% Define the biological optimization model for treatment planning

% if an alpha-beta table is used then define the quantityTable like this
pln.bioModel                                = matRad_bioModel(pln.radiationMode, 'doseAveragedTabulatedAlphaBeta');
pln.bioModel.quantityTableName              = 'AlphaBetaTable_LEMI_Scholz06_AX313_BX615';

% choose the stopping power that is used to average the alpha and beta by
% the dose (stopping power and fluence)
pln.bioModel.stoppingPowerTableName         = 'SPtable';

% select the fragments in the table with their Z and A
selectZ = [1,2,3,4,5,6];
selectA = [1,4,7,9,11,12];

pln.bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);

pln.propOpt.quantityOpt = 'physicalDose';
pln.propDoseCalc.engine = 'HongPB';

%% z* table
% if a z* table is used the table must be called like this
% pln.bioModel = matRad_bioModel(pln.radiationMode, 'doseAveragedTabulatedAlphaBeta');
% pln.bioModel.ZstarTableName = 'Ztable_mMKM_zs_Inaniwa10_A0172_B0615';
% pln.bioModel.stoppingPowerTableName = 'SPtable';
% 
% selectZ = [1,2,3,4,5,6];
% 
% selectA = [1,4,7,9,11,12];
% 
% pln.bioModel.includedFragments = arrayfun(@(z,a) struct('Z', z, 'A', a), selectZ, selectA);
% pln.propOpt.quantityOpt = 'physicalDose';
% pln.propDoseCalc.engine = 'HongPB';

%% Dose calculation
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Inverse Optimization  for Carbon ion treatment based on RBE-weighted dose
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.

resultGUI = matRad_fluenceOptimization(dij,cst,pln);