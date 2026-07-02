%load('carbon_Generic_clusterDose_prestep.mat')
%energies = [machine.data(121).Fluence.energyBin];
%energies(1) = 0.01;
% for i = 1:size(machine.data(121).Fluence.energyBin,2)
%     energies = [energies machine.data(i).energy];
% end

% Topas requires energy per nucleon to be the same for all ions, so fix
% this to be energy per nucleon
%energies = logspace(-1,log10(max(machine.data(121).Fluence.energyBin(:))),250);
energies = logspace(-1,log10(1000),100);


%model = 'LEMIII';
model = 'mMKM';

%calculusType = 'rapidLEM_Scholz2006';
%calculusType = 'rapidLEM_Russo2011';
%calculusType = 'rapidMKM_MCFMKM';
%calculusType = 'rapidMKM_Kase2008';
calculusType = 'rapidmMKM_Inaniwa2010';

survivalParameters.ions = {'H','He', 'Li', 'Be', 'B', 'C','N','O'};
survivalParameters.ionsA = [1,4,7,9,11,12,14,16];

%Define parameters to compute (for the time being fix alpha/beta)
survivalParameters.MKM_alpha0 = [0.172];
survivalParameters.MKM_beta0  = [0.0615];
%survivalParameters.LEM_Dt     = [30];                      %Taken from LEM3 paper, figure 5 (ACCURACY OF THE LOCAL EFFECT MODEL FOR THE PREDICTION OF BIOLOGIC EFFECTS OF CARBON ION BEAMS IN VITRO AND IN VIVO)
%survivalParameters.LEM_rNucleus     = 5*ones(1,size(survivalParameters.LEM_Dt,2));
survivalParameters.MKM_rNucleus = 3.9;
survivalParameters.MKM_rDomain = 0.32;

survivalParameters.survivalParameterFileName = model;
survivalParameters.output          = 'LQ_pars';
survivalParameters.model           = model;
survivalParameters.calculusType    = calculusType;
survivalParameters.parallelismType = '0';
survivalParameters.cellType        = 'HSG';
survivalParameters.trackMode = 'histogram';
survivalParameters.energies  = energies;

survivalParameters.wDir     = fullfile('userdata/survivalTables');

survivalParameters.fileName = 'mMKM_MATLAB';
survivalParameters.survivalSourcePath = '/home/lisa/Survival';


genParameterFile(survivalParameters);

% now manually run the survival code

%% Read out the CSV
%survivalParameters.csvName = [survivalParameters.fileName,'_rapidLEM_Scholz2006_LQparameters_',model,'.csv'];
%survivalParameters.csvName = [survivalParameters.fileName,'_rapidMKM_Kase2008_LQparameters_',model,'.csv'];
survivalParameters.csvName = ['mMKM01_zStar_O_LQparameters_mMKM.csv'];
%[data,meta] = readoutCSV(fullfile('userdata/survivalTables', survivalParameters.csvName));
[data,meta] = readoutCSV(fullfile('userdata/survivalTables',survivalParameters.csvName));

saveDir = fullfile('userdata');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
RBEtable.meta = meta;
RBEtable.data = data;

save(fullfile(saveDir, sprintf('%s_AX%1.2f_BX%1.2f.mat', meta.model, meta.alphaX, meta.betaX)), 'RBEtable');

function genParameterFile(survivalParameters)
         
         if ~exist(survivalParameters.wDir, 'dir')
            mkdir(survivalParameters.wDir);
         end

         fileName = fullfile(survivalParameters.wDir, [survivalParameters.fileName, '.sh']);
         %for n = 1:size(fileName,2)
         fID = fopen(fileName, 'w'); %fopen(['//wsl.localhost/Ubuntu',fileName], 'w','native', 'UTF-8'); % UNC path for writing from Windows
         
         if fID == -1
             error('fopen failed! Could not open file: %s', fileName);
         else
             disp(['Successfully opened file: ', fileName]);
         end
         %disp(['Attempting to save file to: ', fileName]);
         
         % fprintf(fID, '#!/bin/bash\n\n');
         % fprintf(fID,'echo $pwd ');
         % fprintf(fID, '\n');

         %source
         fprintf(fID, 'source ');
         fprintf(fID, survivalParameters.survivalSourcePath);
         fprintf(fID, strcat('/setenv.sh\n'));       
         fprintf(fID, '\n');
         
         %projectName
         fprintf(fID, 'projectName="');
         fprintf(fID, sprintf('%s_%s', survivalParameters.fileName, survivalParameters.calculusType));
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %output
         fprintf(fID, 'output="');
         fprintf(fID, survivalParameters.output);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %model
         fprintf(fID, 'model="');
         fprintf(fID, survivalParameters.model);
         fprintf(fID, '"\n');
         
         %calculusType
         fprintf(fID, 'calculusType="');
         fprintf(fID, survivalParameters.calculusType);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %parallelismType
         fprintf(fID, 'parallelismType="');
         fprintf(fID, survivalParameters.parallelismType);
         fprintf(fID, '" \n');
         fprintf(fID, '\n');
         
         %cellType
         fprintf(fID, 'cellType="');
         fprintf(fID, survivalParameters.cellType);
         fprintf(fID, '" \n');
         fprintf(fID, '\n');
         
         %Model parameters
         switch survivalParameters.model

             case 'MKM'
               fprintf(fID, 'MKM_alpha0=%1.2f', survivalParameters.MKM_alpha0);
               fprintf(fID, '\n');
               
               fprintf(fID, 'MKM_beta0=%1.2f',survivalParameters.MKM_beta0);
               fprintf(fID, '\n');
               
               fprintf(fID, 'MKM_rNucleus=%1.2f', survivalParameters.MKM_rNucleus);
               fprintf(fID, '\n');
               
               fprintf(fID, 'MKM_rDomain=%1.2f',survivalParameters.MKM_rDomain);
               fprintf(fID, '\n');
               fprintf(fID, '\n');
             
             case {'LEMI', 'LEMII', 'LEMIII'}
               fprintf(fID, 'LEM_alpha0=%1.2f', survivalParameters.LEM_alpha0);
               fprintf(fID, '\n');
               
               fprintf(fID, 'LEM_beta0=%1.2f',survivalParameters.LEM_beta0);
               fprintf(fID, '\n');
               
               fprintf(fID, 'LEM_rNucleus=%1.2f', survivalParameters.LEM_rNucleus);
               fprintf(fID, '\n');
               
               fprintf(fID, 'LEM_Dt=%1.2f',survivalParameters.LEM_Dt);
               fprintf(fID, '\n');
               fprintf(fID, '\n');


         end
         
         % Declare array
         fprintf(fID, 'ions=(');
         for i=1:numel(survivalParameters.ions)
             fprintf(fID, '"%s" ', survivalParameters.ions{i});
         end
         fprintf(fID, ')\n\n');

         fprintf(fID, 'ionsA=(');
         for i=1:numel(survivalParameters.ionsA)
             fprintf(fID, '%d ', survivalParameters.ionsA(i));
         end
         fprintf(fID, ')\n\n');
         
         %trackMode
         fprintf(fID, 'trackMode="');
         fprintf(fID,survivalParameters.trackMode);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %energies
         fprintf(fID, 'energies="');
         fprintf(fID, '%3.4f ', survivalParameters.energies);
         fprintf(fID, '"\n');
         
         fprintf(fID, '\n');

         fprintf(fID, 'for ((i = 0; i < ${#ions[@]}; i++)); do\n');
         fprintf(fID, '\t ion=${ions[i]}\n');
         fprintf(fID, '\t ionA=${ionsA[i]}\n');
         fprintf(fID, '\t echo "Running for ion: $ion"\n\n');
         fprintf(fID, '\t # Multiply all energies at once using awk and capture output\n');
         fprintf(fID, '\t total_energies=$(echo "$energies" | awk -v a="$ionA" ''{\n');
         fprintf(fID, '\t \tfor (i = 1; i <= NF; i++) {\n');
         fprintf(fID, '\t \t \t printf "%%.6f ", $i * a;\n');
         fprintf(fID, '\t \t}\n');
         fprintf(fID, '\t }'')\n\n');
         fprintf(fID, '\t # Remove trailing space\n');
         fprintf(fID, '\t total_energies=${total_energies%%%% }\n');

         
         fprintf(fID, '\n');
         
         fprintf(fID, [survivalParameters.survivalSourcePath, '/CMakeProject     -projectName $projectName\\\n']);
         fprintf(fID, '                -output $output \\\n');
         fprintf(fID, '                -model $model \\\n');
         fprintf(fID, '                -calculusType $calculusType \\\n');
         fprintf(fID, '                -parallelismType $parallelismType \\\n');
         fprintf(fID, '                -cellType $cellType \\\n');
         %fprintf(fID, '                -LEM_alpha0 $LEM_alpha0 \\\n');
         %fprintf(fID, '                -LEM_beta0 $LEM_beta0 \\\n');
         %fprintf(fID, '                -LEM_rNucleus $LEM_rNucleus \\\n');
         %fprintf(fID, '                -LEM_Dt $LEM_Dt \\\n');
         fprintf(fID, '                -MKM_alpha0 $MKM_alpha0 \\\n');
         fprintf(fID, '                -MKM_beta0 $MKM_beta0 \\\n');
         fprintf(fID, '                -MKM_rNucleus $MKM_rNucleus \\\n');
         fprintf(fID, '                -MKM_rDomain $MKM_rDomain \\\n');
         fprintf(fID, '                -ion $ion \\\n');
         fprintf(fID, '                -trackMode $trackMode \\\n');
         fprintf(fID, '                -energies $energies \\;\n');
         fprintf(fID, 'done\n');
         
         fclose(fID);
      
end

function [data, meta] = readoutCSV(filename)

    rawData = readtable(filename);

    data.alphaX = rawData.alpha_0(1);
    data.betaX  = rawData.beta_0(1);

    ionsNames = rawData.particle;
    uniqueIonNames = unique(ionsNames, 'stable');

    %data.energies(:,i) = rawData.meanEnergy(strcmp(ionsNames, uniqueIonNames{1}));
    % data.energies = data.energies(:)'; % Force row array

    for i=1:numel(uniqueIonNames)
        ionIdx = strcmp(ionsNames, uniqueIonNames{i});
        data.alpha(:,i) = rawData.alpha(ionIdx);
        data.beta(:,i)  = rawData.beta(ionIdx);
        %data.zs(:,i) = rawData.z(ionIdx);
        data.energies(:,i) = (rawData.meanEnergy(ionIdx))';

        switch uniqueIonNames{i}
            case 'H'
                data.includedIons(i).Z = 1;
                data.includedIons(i).A = 1;
            case 'He'
                data.includedIons(i).Z = 2;
                data.includedIons(i).A = 4;
            case 'Li'
                data.includedIons(i).Z = 3;
                data.includedIons(i).A = 7;
            case 'Be'
                data.includedIons(i).Z = 4;
                data.includedIons(i).A = 9;
            case 'B'
                data.includedIons(i).Z = 5;
                data.includedIons(i).A = 11;
            case 'C'
                data.includedIons(i).Z = 6;
                data.includedIons(i).A = 12;
            case 'N'
                data.includedIons(i).Z = 7;
                data.includedIons(i).A = 14;
            case 'O'
                data.includedIons(i).Z = 8;
                data.includedIons(i).A = 16;
        end        
    end

    meta.model           = sprintf('%s_%s', rawData.model{1}, rawData.calculusType{1});
    meta.description     = 'Model created with survival'; 
    meta.modelParameters.alphaX =   data.alphaX;
    meta.modelParameters.betaX  =   data.betaX;
    meta.modelParameters.rNucleus = rawData.r_nucleus(1); 
    %meta.modelParameters.Dt     = rawData.D_t(1);
    meta.modelParameters.rDomain = rawData.r_domain(1);
    meta.alphaX          = data.alphaX;
    meta.betaX           = data.betaX;


end

%% Read out csv for dE/dx
load ("RBEtable_rapidLEM_Scholz2006_allIons_LEMIII0102_new.mat");
dEdx.csvName = ['sp_table_water_icru.csv'];
for i = 1:10
    RBEtable = readout_CSV(fullfile('userdata/survivalTables', dEdx.csvName),i);
end
% for i = 1:6
%     switch RBEtable.data(1).includedIons(i).Z
%         case 1
%             dEdx.csvName = ['dE-dx_H_logEnergy_crange.csv'];
%             RBEtable = readout_CSV(fullfile('userdata/survivalTables', dEdx.csvName),i);
%         case 2
%             dEdx.csvName = ['dE-dx_He_logEnergy_crange.csv'];
%             RBEtable = readout_CSV(fullfile('userdata/survivalTables', dEdx.csvName),i);
%         case 3
%             dEdx.csvName = ['dE-dx_Li_logEnergy_crange.csv'];
%             RBEtable = readout_CSV(fullfile('userdata/survivalTables', dEdx.csvName),i);
%         case 4
%             dEdx.csvName = ['dE-dx_Be_logEnergy_crange.csv'];
%             RBEtable = readout_CSV(fullfile('userdata/survivalTables', dEdx.csvName),i);
%         case 5
%             dEdx.csvName = ['dE-dx_B_logEnergy_crange.csv'];
%             RBEtable = readout_CSV(fullfile('userdata/survivalTables', dEdx.csvName),i);
%         case 6
%             dEdx.csvName = ['dE-dx_C_logEnergy_crange.csv'];
%             RBEtable = readout_CSV(fullfile('userdata/survivalTables', dEdx.csvName),i);
%     end
% end

function RBEtable = readout_CSV(filename,i)

    rawData = readtable(filename);
    
    load ("RBEtable_rapidLEM_Scholz2006_allIons_LEMIII0102_new.mat");
    
    ionName = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'};

    RBEtable.data(1).dEdx(:,i) = rawData.(ionName{i});
    RBEtable.data(2).dEdx(:,i) = rawData.(ionName{i});

    RBEtable.meta.dEdxCalculator     = 'https://aptg.github.io/libdedx/'; 
    
    save("RBEtable_rapidLEM_Scholz2006_allIons_LEMIII0102_new.mat","RBEtable");

end
% function RBEtable = readout_CSV(filename,i)
% 
%     rawData = readtable(filename);
% 
%     load ("RBEtable_rapidMKM_Kase2008_allIons_MKM_0102.mat");
% 
%     RBEtable.data(1).dEdx(:,i) = rawData.Result;
%     RBEtable.data(2).dEdx(:,i) = rawData.Result;
% 
%     RBEtable.meta.dEdxCalculator     = 'dE/dx calculated with B.A. Weaver Crange package'; 
% 
%     save("RBEtable_rapidMKM_Kase2008_allIons_MKM_0102.mat","RBEtable");
% 
% end

%% Read out the energy of dEdx
dEdx.csvName = ['sp_table_water_icru.csv'];
dataCSV = readtable(fullfile('userdata/survivalTables',dEdx.csvName));

for i = 1:10
    SPtable.data(i).energies = dataCSV.EnergyPerNucleon_MeV_nucleon_;
end

%%
%NEW

tableMCF = readtable('userdata/survivalTables/MCF_MKM_Shannon_LQparameters_MKM.csv');
tableKase = readtable('userdata/survivalTables/KaseMKM_Shannon_LQparameters_MKM.csv');

%% 
figure

AlphaC_MCF = tableMCF.alpha(end-99:end);
meanEnergyC_MCF = tableMCF.meanEnergy(end-99:end);
meanEnergyC_MCF_MeVu = meanEnergyC_MCF / 12;
AlphaC_Kase = tableKase.alpha(end-97:end);
meanEnergyC_Kase = tableKase.meanEnergy(end-97:end);
meanEnergyC_Kase_MeVu = meanEnergyC_Kase / 12;

loglog(meanEnergyC_MCF_MeVu, AlphaC_MCF, '.-', 'DisplayName', 'MCF MKM');

% figure;
hold on;
% loglog(RBEtable.data(6).energies, RBEtable.data(6).alpha(:,1), '.-');
loglog(meanEnergyC_Kase_MeVu, AlphaC_Kase, '.-', 'DisplayName','MKM Kase2008');

legend();
grid on;

ylabel('\alpha [Gy^{-1}]');
xlabel('Energy [MeV/u]');


figure

BetaC_MCF = tableMCF.beta(end-99:end);
BetaC_Kase = tableKase.beta(end-97:end);

loglog(meanEnergyC_MCF_MeVu, BetaC_MCF, '.-', 'DisplayName', 'MCF MKM');

% figure;
hold on;
% loglog(RBEtable.data(6).energies, RBEtable.data(6).alpha(:,1), '.-');
loglog(meanEnergyC_Kase_MeVu, BetaC_Kase, '.-', 'DisplayName','MKM Kase2008');

legend();
grid on;

ylabel('\beta [Gy^{-2}]');
xlabel('Energy [MeV/u]');

%%
figure
loglog(tableMCF.meanEnergy(:), tableMCF.beta(:), '.-', 'DisplayName', 'MCF MKM');
hold on
loglog(tableKase.meanEnergy(:), tableKase.beta(:), '.-', 'DisplayName','MKM Kase2008');

legend();
grid on;

ylabel('\beta [Gy^{-2}]');
xlabel('Energy [MeV/u]');
