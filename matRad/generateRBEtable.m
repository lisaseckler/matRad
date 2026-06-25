function RBEtableName = generateRBEtable(pln,model,calculusType,alpha0,beta0,Dt,RN,SurvivalPath,folder,survivalCommand,FilePath)
%Generate RBEtables interfacing with survival code
matRad_cfg = MatRad_Config.instance();
%Define energies (MeV/u)
%energies = logspace(-1,log10(510),250);
load(['basedata\',pln.radiationMode,'_',pln.machine,'.mat']);
energies = [];
for i = 1:size(machine.data,2)
    energies = [energies machine.data(i).energy];
end


%Define ions
% ions  = {'H', 'He', 'Li', 'Be', 'B', 'C'};
% ionsMassNumber = [1 4 7 9 11 12];
ions  = {'H', 'He', 'Li', 'Be', 'B','C','N','O'};
ionsMassNumber = [1.0079,4.0026,6.941,9.0122,10.811,12.0107,14.0067,15.9994];

%Define model and calculus type
if strcmp(model, 'MKM')
   %model = 'MKM'; % MKM is not working right now!
   base = 'MKM';
   modelSpecific = 'rDomain';
   rDomain = Dt;
   rNucleus = RN;
else
   base = 'LEM';
   modelSpecific = 'Dt';
end

%calculusType = 'rapidLEM_Scholz2006';
%calculusType = 'rapidMKM_Kase2008';

%Define parameters to compute (for the time being fix alpha/beta)

% alpha0 = [0.10];
% beta0  = [0.05];

%model specific

%Dt     = [30];                      %Taken from LEM3 paper, figure 5 (ACCURACY OF THE LOCAL EFFECT MODEL FOR THE PREDICTION OF BIOLOGIC EFFECTS OF CARBON ION BEAMS IN VITRO AND IN VIVO)
%RN     = 8*ones(1,size(Dt,2));
% MKM_alpha0=0.10
% MKM_beta0=0.05
% rNucleus = 3.90;
% rDomain = 0.38;

%Define the interface class
IntegralS = integralDoseSurvival();
IntegralS.wDir = [fullfile(SurvivalPath,folder)];

%This is source path see from wsl
IntegralS.survivalSourcePath = SurvivalPath;
IntegralS.survivalCommand = survivalCommand;

IntegralS.survivalParameterFileName = model;
IntegralS.calcProperties.output          = 'LQ_pars';
IntegralS.calcProperties.model           = model;
IntegralS.calcProperties.calculusType    = calculusType;
%IntegralS.calcProperties.precision       = 0.15;
IntegralS.calcProperties.parallelismType = '0';
IntegralS.calcProperties.cellType        = 'HSG';

IntegralS.calcProperties.trackMode = 'histogram';
IntegralS.calcProperties.energies  = energies;
IntegralS.calcProperties.ion = {};

% idx = strfind(IntegralS.wDir, '\');
% dirPath = IntegralS.wDir(idx(1)+1:end);
% dirPath(strfind(dirPath, '\')) = '/';
% IntegralS.wDirWsl = ['/home/lisa/CMakeProject/'];

%Model parameters
switch base

    case 'MKM'
        for alphaIdx = 1:size(alpha0,2)
            IntegralS.calcProperties.modelParam.MKM_alpha0   = alpha0(alphaIdx);
            IntegralS.calcProperties.modelParam.MKM_beta0    = beta0(alphaIdx);

            for pIdx=1:size(rDomain,2)

                IntegralS.calcProperties.projectName     = ['LUT_alpha_',num2str(alpha0(alphaIdx)),'_rDomain_', num2str(rDomain(pIdx))];

                IntegralS.calcProperties.modelParam.MKM_rNucleus = rNucleus(pIdx);
                IntegralS.calcProperties.modelParam.MKM_rDomain  = rDomain(pIdx);
                matRad_cfg.dispInfo('alpha value %d, rDomain value %d \n', alpha0(alphaIdx), rDomain(pIdx));
                for k=1:size(ions,2)
                    IntegralS.calcProperties.energies  = energies.*ionsMassNumber(k);
                    IntegralS.calcProperties.ion{k} = ions{k};
                end
                IntegralS.genParameterFile();
                filePath = strrep(strcat(IntegralS.wDir, filesep,IntegralS.survivalParameterFileName, '.txt'),'\','/');
                %IntegralS.survivalExecutionCommand();
                %IntegralS.savingCommand();
                % Run the script inside WSL
                IntegralS.survivalExecutionCommand = ['wsl bash ',filePath];

                %Execution with system does not work, check -> Solved, if does not work
                %again, check that default destribution is Ubuntu and not Ubunutu-20.04.
                %Else, from prompt wsl --setdefault Ubuntu
                matRad_cfg.dispInfo('Executing Survival... this may take a while...');
                a = IntegralS.execute();
                if a == 0
                    matRad_cfg.dispInfo('done. \n');
                else
                    matRad_cfg.dispInfo('error. \n');
                end
            end
        end
        case 'LEM'
        for alphaIdx = 1:size(alpha0,2)
            IntegralS.calcProperties.modelParam.LEM_alpha0   = alpha0(alphaIdx);
            IntegralS.calcProperties.modelParam.LEM_beta0    = beta0(alphaIdx);

            for pIdx=1:size(Dt,2)

                IntegralS.calcProperties.projectName     = ['LUT_alpha_',num2str(alpha0(alphaIdx)),'_Dt_', num2str(Dt(pIdx))];

                IntegralS.calcProperties.modelParam.LEM_rNucleus = RN(pIdx);
                IntegralS.calcProperties.modelParam.LEM_Dt  = Dt(pIdx);
                matRad_cfg.dispInfo('alpha value %d, Dt value %d \n', alpha0(alphaIdx), Dt(pIdx));
                for k=1:size(ions,2)
                    IntegralS.calcProperties.energies  = energies.*ionsMassNumber(k);
                    IntegralS.calcProperties.ion{k} = ions{k};
                end
                IntegralS.genParameterFile();
                filePath = strrep(strcat(IntegralS.wDir, filesep,IntegralS.survivalParameterFileName, '.txt'),'\','/');
                %IntegralS.survivalExecutionCommand();
                %IntegralS.savingCommand();
                % Run the script inside WSL
                IntegralS.survivalExecutionCommand = ['wsl bash ',filePath];

                %Execution with system does not work, check -> Solved, if does not work
                %again, check that default destribution is Ubuntu and not Ubunutu-20.04.
                %Else, from prompt wsl --setdefault Ubuntu
                matRad_cfg.dispInfo('Executing Survival... this may take a while...');
                a = IntegralS.execute();
                if a == 0
                    matRad_cfg.dispInfo('done. \n');
                else
                    matRad_cfg.dispInfo('error. \n');
                end
            end
        end
end

%% read Out
clear RBEtable;
%meta info
switch base

    case 'MKM'
        filenames = cell(1,numel(model));
        for pIdx = 1:size(rDomain,2)

            for n = 1:size(model,2)
                entryName = sprintf('entry%d', n);
            RBEtable.(entryName).meta.model = [calculusType, '_',modelSpecific,'_', num2str(rDomain(pIdx))];
            RBEtable.(entryName).meta.description = ['Model obtained with', calculusType, ' computed with Survival code.'];

            modelParameters = struct('alphaX', [], ...
                             'betaX', [], ...
                             'rNucleus', [], ...
                             'Dt', []);

            for alphaIdx =1:size(alpha0,2)
                modelParameters.alphaX = alpha0;
                modelParameters.betaX  = beta0;
                modelParameters.rNucleus = rNucleus(pIdx);
                modelParameters.rDomain  = rDomain(pIdx);

                RBEtable.(entryName).meta.modelParameters = modelParameters;

                IntegralS.exelFile = [fullfile(FilePath,'matRad')];
                %filenames = dir([IntegralS.wDir, filesep, '*.csv']);
                for l = 1:size(model,2)
                    filenames{l} = {[strrep([IntegralS.exelFile, filesep,'LUT_alpha_',num2str(alpha0(alphaIdx)),'_',modelSpecific,'_', num2str(Dt(pIdx)), '_LQparameters_',model{l},'.csv'],'\','/')]};
                end
                [alphaE, betaE] = IntegralS.readMultipleIonLUT(ions, filenames);

                %This is directly alpha, not zD;
                %         for iIdx = 1:size(ions)
                %             varAlpha(:,k) = alpha0(alphaIdx) + beta0(aphaIdx) * alphaE;
                %             varBeta(:,k) = betaE;
                %         end

                RBEtable.(entryName).data(alphaIdx).alphaX = alpha0(alphaIdx);
                RBEtable.(entryName).data(alphaIdx).betaX = beta0(alphaIdx);


                RBEtable.(entryName).data(alphaIdx).energies     = energies;

                RBEtable.(entryName).data(alphaIdx).includedIonZ = [1 2 3 4 5 6 7 8];

                RBEtable.(entryName).data(alphaIdx).alpha = struct();  % Initialize struct
                for m = 1:size(alphaE,2)
                    ionName = ions{m};
                    RBEtable.(entryName).data(alphaIdx).alpha.(ionName) = alphaE{m};
                    RBEtable.(entyName).data(alphaIdx).beta.(ionName) = betaE{m};
                end

                %         RBEtable.data(alphaIdx).alpha        = alphaE;
                %         RBEtable.data(alphaIdx).beta         = betaE;
            end

            RBEtableName = [matRad_cfg.matRadRoot, filesep, 'matRad', filesep, 'bioModels', filesep, 'RBEtables', filesep,'RBEtable_', calculusType, '_allIons_',model{n},'_', num2str(rDomain(pIdx)), '.mat'];
            save(RBEtableName, 'RBEtable');
            
            end
            RBEtable = [];
        end

        case 'LEM'
            filenames = cell(1,numel(model));
        for pIdx = 1:size(Dt,2)
            % for both not only for one... still not done
            for n = 1:size(model,2)
                entryName = sprintf('entry%d', n);
            RBEtable.(entryName).meta.model = [calculusType, '_',modelSpecific,'_', num2str(Dt(pIdx))];
            RBEtable.(entryName).meta.description = ['Model obtained with', calculusType, ' computed with Survival code.'];

            modelParameters = struct('alphaX', [], ...
                             'betaX', [], ...
                             'rNucleus', [], ...
                             'Dt', []);

            for alphaIdx =1:size(alpha0,2)
                modelParameters.alphaX = alpha0;
                modelParameters.betaX  = beta0;
                modelParameters.rNucleus = RN(pIdx);
                modelParameters.Dt  = Dt(pIdx);

                RBEtable.(entryName).meta.modelParameters = modelParameters;

                IntegralS.exelFile = [fullfile(FilePath,'matRad')];
                %filenames = dir([IntegralS.wDir, filesep, '*.csv']);
                
                for l = 1:size(model,2)
                    filenames{l} = {[strrep([IntegralS.exelFile, filesep,'LUT_alpha_',num2str(alpha0(alphaIdx)),'_',modelSpecific,'_', num2str(Dt(pIdx)), '_LQparameters_',model{l},'.csv'],'\','/')]};
                end
                [alphaE, betaE] = IntegralS.readMultipleIonLUT(ions, filenames);

                %This is directly alpha, not zD;
                %         for iIdx = 1:size(ions)
                %             varAlpha(:,k) = alpha0(alphaIdx) + beta0(aphaIdx) * alphaE;
                %             varBeta(:,k) = betaE;
                %         end

                RBEtable.(entryName).data(alphaIdx).alphaX = alpha0(alphaIdx);
                RBEtable.(entryName).data(alphaIdx).betaX = beta0(alphaIdx);


                RBEtable.(entryName).data(alphaIdx).energies     = energies;

                RBEtable.(entryName).data(alphaIdx).includedIonZ = [1 2 3 4 5 6 7 8];

                RBEtable.(entryName).data(alphaIdx).alpha = struct();  % Initialize struct
                for m = 1:size(alphaE,2)
                    ionName = ions{m};
                    RBEtable.(entryName).data(alphaIdx).alpha.(ionName) = alphaE{m};
                    RBEtable.(entryName).data(alphaIdx).beta.(ionName) = betaE{m};
                end

                %         RBEtable.data(alphaIdx).alpha        = alphaE;
                %         RBEtable.data(alphaIdx).beta         = betaE;
            
            data = RBEtable.(entryName);
            RBEtableName = [matRad_cfg.matRadRoot, filesep, 'matRad', filesep, 'bioModels', filesep, 'RBEtables', filesep,'RBEtable_', calculusType, '_allIons_',model{n},'_', num2str(Dt(pIdx)), [entryName, '.mat']];
            save(RBEtableName, 'data');
            
            end
            end
            RBEtable = [];
        end
end

%% Print to topas table (To Be Tested!)
for n = 1:size(model,2)
    entryName = sprintf('entry%d', n);
switch base
    case 'MKM'
        RBEtable.(entryName) = load([FilePath,'\matRad\bioModels\RBEtables\RBEtable_',calculusType,'_allIons_',IntegralS.survivalParameterFileName{n},'_',num2str(IntegralS.calcProperties.modelParam.MKM_rDomain),'.mat']);
    case 'LEM'
        RBEtable.(entryName) = load([FilePath,'\matRad\bioModels\RBEtables\RBEtable_',calculusType,'_allIons_',IntegralS.survivalParameterFileName{n},'_',num2str(IntegralS.calcProperties.modelParam.LEM_Dt),'.mat']);
end
end
for n = 1:size(model,2)
    entryName = sprintf('entry%d', n);
    fID = fopen(['carbonTable_',model{n},'.txt'], 'w');
    fprintf(fID, '### Survival generated table ###\n');
    fprintf(fID, 'sv:Sc/CellGeneric_abR2/HCP/ParticleName 		= 6 "Proton" "Helium" "Lithium" "Beryllium" "Boron" "Carbon"\n');
    fprintf(fID, 'iv:Sc/CellGeneric_abR2/HCP/ParticleZ    		= 6 1 2 3 4 5 6\n');
    fprintf(fID, 'dv:Sc/CellGeneric_abR2/HCP/KineticEnergyPerNucleon 	= %i',length(RBEtable.(entryName).data.energies));
    fprintf(fID, ' %1.3f', RBEtable.(entryName).data.energies);
    fprintf(fID, ' MeV \n');
    fprintf(fID, '\n');

    ionsString = {'Proton', 'Helium', 'Lithium', 'Beryllium', 'Boron', 'Carbon'};
    for k=1:6
        ionName = ions{k};
        fprintf(fID, 'dv:Sc/CellGeneric_abR2/HCP/%s/Alpha 	= %i',ionsString{k},length(RBEtable.(entryName).data.energies));
        fprintf(fID, ' %1.3f', RBEtable.(entryName).data.alpha.(ionName));
        fprintf(fID, ' /Gy\n');
    end
    fprintf(fID, '\n');
    for k=1:6
        ionName = ions{k};
        fprintf(fID, 'dv:Sc/CellGeneric_abR2/HCP/%s/Beta 	= %i',ionsString{k},length(RBEtable.(entryName).data.energies));
        fprintf(fID, ' %1.3f', RBEtable.(entryName).data.beta.(ionName));
        fprintf(fID, ' /Gy2\n');
    end
    fclose(fID);
end
end