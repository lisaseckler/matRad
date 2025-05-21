matRad_rc;
matRad_cfg = MatRad_Config.instance();

load('BOXPHANTOM.mat');

pln.radiationMode = 'carbon';
pln.machine       = 'Generic_clusterDose_preStep_mod';
pln.propDoseCalc.calcLET = 0;

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);

pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 8;
pln.propDoseCalc.doseGrid.resolution.y = 8;
pln.propDoseCalc.doseGrid.resolution.z = 8;

pln.multScen = matRad_multScen(ct, 'nomScen');

%% stf
stf = matRad_generateStf(ct,cst,pln);

%% Set cst
cst{1,6} = {};
cst{2,6} = {};

cst{1,6} = {struct(DoseObjectives.matRad_SquaredOverdosing(5,4.2))};
cst{1,6}{1}.robustness = 'none';
cst{1,6}{1}.quantity   = 'physicalDose';

cst{2,6} = {struct(DoseObjectives.matRad_SquaredDeviation(100,42))};
cst{2,6}{1}.robustness = 'none';
cst{2,6}{1}.quantity   = 'physicalDose';

%% Dose calc

% Select biological model
% % Available models are:
% %   RBEminMax (LET based): MCN, WED, CAR, LSM, (protons)
% %                          HEL                 (helium)
% %   kernel based:          LEM                 (carbon)
% %   Tabulated              TAB                  (proton, helium, carbon) (requires RBEtable and spectra in base data)
% 
% % Example 1: load the MCN model. The machine input is optional in this case.
% % The output is an instance of the model. 
% pln.bioModel = matRad_bioModel(pln.radiationMode,'MCN', 'generic_MCsquare');
% 
% % Alternatively, the biological model can also be assigned as a structure.
% % The only andatory field in this case is the 'model' field:
% % pln.bioModel.model = 'MCN';
% 
% % Example 2: Some models have parameters that can be tuned by the user.
% % For example, we can instantiate a constRBE model
% %
% % pln.bioModel = matRad_bioModel(pln.radiationMode, 'constRBE');
% %
% % and assign a custom value for the constant RBE
% %
% % pln.bioModel.RBE = 4.2;


% Example 3: Tabulated models require  the additional specification of an
% RBE table and additional base data entries. Additional details are
% provided in the Wiki
pln.bioModel = matRad_bioModel(pln.radiationMode, 'TAB');
pln.bioModel.fragmentsToInclude  = {'H1','C'};
pln.bioModel.RBEtableName       = 'RBEtable_rapidLEM_Scholz2006_allIons_LEMI_30';

% ionName = fields(RBEtable.data.alpha);
% RBEtable.data.alpha1 = [];
% for i = 1:size(ionName,1)
%     RBEtable.data.alpha1 = [RBEtable.data.alpha1 RBEtable.data.alpha.(ionName{i})];
% end
 
pln.propOpt.quantityOpt = 'RBExDose';
pln.propDoseCalc.engine = 'HongPB';
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Fluence opt

resultGUI = matRad_fluenceOptimization(dij,cst,pln);


%matRadGUI;

%% asd
figure
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI.effect,3,80);
figure
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI.physicalDose,3,80);
figure
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_Topas.effect,3,80);
figure
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_Topas.physicalDose,3,80);


%% 
pln.propDoseCalc.engine = 'TOPAS';
pln.propDoseCalc.externalCalculation = 'write';
pln.propDoseCalc.scorer.calcBioDose = true;


resultGUI_Topas = matRad_calcDoseDirect(ct,stf,pln,cst, resultGUI.w);

