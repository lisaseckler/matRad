load carbon_Generic_clusterDose_prestep_update.mat;
load('RBEtable_rapidMKM_Kase2008_allIons_MKM_0102.mat');

%% Load zStar
load DataInaniwa2010.mat
mMKM_zValues = readtable('\\wsl.localhost\Ubuntu\home\niklas\dev\Survival\mMKM_zs_mMKM.csv');
Zs = 1:6;
strings = ["H","He","Li","Be","B","C"];
As = [1 4 7 9 10 12];

dedx_icru = readtable('\\wsl.localhost\Ubuntu\home\niklas\dev\Survival\sp_table_water_icru.csv');

%%
particle = [1 4 5 6 7 8];
fields = fieldnames(RBEtable);
for i = 1:numel(fields)
    tableEntry.data(i) = RBEtable.(fields{i,1});
end

dEdx = SPtable.data(6).dEdx;
%modelEnergies = tableEntry.data(1).data(1).energies;
modelEnergies = RBEtable.data(6).energies; % MeV

bdEntry = machine.data(48);
% bdEntry = LEM.data(61);


depths = bdEntry.depths;
depthGrad = gradient(depths);
fluence = machine.data(48).Fluence;

energies = fluence.spectra(8).energyBin;

%%
figure;
hDeDx = subplot(2,2,1);
hZstar = subplot(2,2,2);
dividend_z = zeros(size(depths));
denominator_z = zeros(size(depths));
dividend_alpha = zeros(size(depths));
denominator_alpha = zeros(size(depths));
dividend_sqrtbeta = zeros(size(depths));
denominator_sqrtbeta = zeros(size(depths));
dividend_LET = zeros(size(depths));
denominator_LET = zeros(size(depths));

mMKM_zValues = Ztable.data;

dEdxEnergy = SPtable.data(6).energies;
modelEnergiesMev = modelEnergies/u;
idx = find(strcmp(fields, 'LEMIII014'));

    %get z star
    spectrum = fluence.spectra;
    currZ = spectrum(particle(6)).Z;
    if currZ > 0 && currZ <= numel(As)
        rowIndices = ismember([mMKM_zValues.Z], currZ);
    else
        rowIndices = [];
    end

    if ~isempty(rowIndices)
        u = As(currZ);
        zEnergies = mMKM_zValues(rowIndices).energies/u; %MeV/u
        zValues = mMKM_zValues(rowIndices).zs;

        specEnergies = spectrum(particle(6)).energyBin; %MeV/u

        dEdxInterp_rbetable = interp1(dEdxEnergy,dEdx(:),specEnergies,"linear","extrap");
        %dEdxInterp = interp1(SPtable(currZ).dEdx,dedx_icru(currZ+1).dEdx,specEnergies,"linear","extrap");
        zInterp = interp1(zEnergies,zValues,specEnergies,"linear","extrap");

        % alphaInterp = interp1(modelEnergiesMev,tableEntry.data(idx).data(currZ).alpha(:,1),specEnergies,"linear","extrap");
        alphaInterpIII014 = interp1(modelEnergies,RBEtable.(fields{idx}).data(6).alpha(:,1),specEnergies,"linear","extrap");
        betaInterpIII014  = interp1(modelEnergies,tableEntry.data(idx).data(currZ).beta(:,1),specEnergies,"linear","extrap");

        hC = semilogx(hDeDx,specEnergies,dEdxInterp,'LineWidth',2.5,'DisplayName',sprintf('Z = %d, A = %d',currZ,As(currZ))); hold(hDeDx,"on");
        semilogx(hDeDx,specEnergies,dEdxInterp_rbetable,'LineWidth',2.5,'DisplayName',sprintf('Z = %d, A = %d',currZ,As(currZ)),'LineStyle',':','Color',hC.Color); hold(hDeDx,"on");
        semilogx(hZstar,specEnergies,zInterp,'LineWidth',2.5,'DisplayName',sprintf('Z = %d, A = %d',currZ,As(currZ))); hold(hZstar,"on");

        %dEdxInterp = dEdxInterp_rbetable;


        zDoseWeighted = zInterp.*dEdxInterp_rbetable.*spectrum(particle(6)).fluenceSpectrum';

        dividend_z = dividend_z + sum(zInterp.*dEdxInterp_rbetable.*spectrum(particle(6)).fluenceSpectrum',2);
        denominator_z = denominator_z + sum(dEdxInterp_rbetable.*spectrum(particle(6)).fluenceSpectrum',2);
        
        dividend_LET = dividend_LET + sum(dEdxInterp_rbetable.^2.*spectrum(particle(6)).fluenceSpectrum',2);
        denominator_LET = denominator_LET + sum(dEdxInterp_rbetable.*spectrum(particle(6)).fluenceSpectrum',2);
        
        dividend_alpha = dividend_alpha + sum(alphaInterpIII014.*dEdxInterp_rbetable.*spectrum(particle(6)).fluenceSpectrum',2);
        denominator_alpha = denominator_alpha + sum(dEdxInterp_rbetable.*spectrum(particle(6)).fluenceSpectrum',2);

        dividend_sqrtbeta = dividend_sqrtbeta + sum(sqrt(betaInterp).*dEdxInterp_rbetable.*spectrum(particle(6)).fluenceSpectrum',2);
        denominator_sqrtbeta = denominator_sqrtbeta + sum(dEdxInterp_rbetable.*spectrum(particle(6)).fluenceSpectrum',2);

    end
%%
legend(hDeDx);
legend(hZstar);

ylabel(hDeDx,'LET / S');
xlabel(hDeDx,'energy [MeV/u]');
grid(hDeDx,"minor");
title(hDeDx,sprintf('Stopping Power (E_{prim} = %g MeV/u)',bdEntry.energy),'interpolated on spectrum energies');

ylabel(hZstar,'z^*');
xlabel(hZstar,'energy [MeV/u]');
grid(hZstar,"minor");
title(hZstar,sprintf('z^* (E_{prim} = %g MeV/u)',bdEntry.energy),'interpolated on spectrum energies');

hMix = subplot(2,2,3);

zMix = dividend_z ./ denominator_z;
zMix(~isfinite(zMix)) = 0;
yyaxis(hMix,"left");
plot(hMix,depths,zMix,'LineWidth',2.5); hold on;
if exist('zDepthC','var')
    plot(hMix,zDepthC.depth,zDepthC.z,'LineWidth',1,'LineStyle',':','Color','k');
end
ylim(hMix,[0 22]);
ylabel(hMix,'z^*');
yyaxis(hMix,"right");
plot(hMix,depths,bdEntry.Z,'LineWidth',1.5,'LineStyle',':');
ylabel(hMix,'IDD');
xlim(hMix,[0 200]);
xlabel(hMix,'depth [mm]');
grid(hMix,"minor");
title(hMix,sprintf('Averaged z* (E_{prim} = %g MeV/u)',bdEntry.energy));

hMixLET = subplot(2,2,4);
LET = dividend_LET ./ denominator_LET;
LET(~isfinite(LET)) = 0;
yyaxis(hMixLET,"left");
plot(hMixLET,depths,LET,'LineWidth',2.5);
ylim(hMixLET,[0 300]);
ylabel(hMixLET,'LET');
yyaxis(hMixLET,"right");
plot(hMixLET,depths,bdEntry.Z,'LineWidth',1.5,'LineStyle',':');
xlim(hMixLET,[0 200]);
ylabel(hMixLET,'IDD');
xlabel(hMixLET,'depth [mm]');
grid(hMixLET,"minor");
title(hMixLET,sprintf('Averaged LET (E_{prim} = %g MeV/u)',bdEntry.energy));

%% Plot alpha
figure;
%alpha0 = 0.0708; %0.282;
%beta0 = 0.0615;
alpha0 = tableEntry.alphaX;
beta0  = tableEntry.betaX;

alphaMKM = (alpha0 + beta0*zMix);
alphaSpec = dividend_alpha ./ denominator_alpha;

plot(depths,alphaMKM,'LineWidth',2.5); hold on;
plot(depths,alphaSpec,'LineWidth',2.5,'LineStyle',':');

%% only one particle

% fields = fieldnames(RBEtable);
%machineDepths = machine.data(48).depths;

dividend_alpha_zero = zeros(size(depths));
denominator_alpha_zero = zeros(size(depths));
dividend_sqrtbeta_zero = zeros(size(depths));
denominator_sqrtbeta_zero = zeros(size(depths));

MCFdepths = MCF_machine.data(76).depths;
MCFalpha = MCF_machine.data(76).alpha;
MCFbeta = MCF_machine.data(76).beta;

% specEnergies = cell(1,6);
% dEdxEnergy = cell(1,6);
% dEdx = cell(1,6);
% dEdxInterp_rbetable = cell(1,6);
% alphaInterp = cell(1,6);
% betaInterp = cell(1,6);
% dividend_alpha_table = cell(1,6);
% dividend_sqrtbeta_table = cell(1,6);
% alphaSpec = cell(1,6);
% sqrtbetaSpec = cell(1,6);

%for p = 1:6
    spectrum = fluence.spectra;
    currZ = spectrum(particle(6)).Z;
    u = As(currZ);
    dEdxEnergy_single = SPtable.data(6).energies;
    dEdx_single = SPtable.data(6).dEdx;
    specEnergies_single = spectrum(particle(6)).energyBin; %MeV/u
    dEdxInterp_rbetable_single = interp1(dEdxEnergy_single,dEdx_single,specEnergies_single,"linear","extrap");

    denominator_alpha_table_single = denominator_alpha_zero + sum(dEdxInterp_rbetable_single.*spectrum(particle(6)).fluenceSpectrum',2); % sum over energy
    %denominator_alpha_table_single = denominator_alpha_zero + sum(spectrum(particle(6)).fluenceSpectrum',2); % sum over energy
    denominator_sqrtbeta_table_single = denominator_sqrtbeta_zero + sum(dEdxInterp_rbetable_single.*spectrum(particle(6)).fluenceSpectrum',2);

    %for i = 1:numel(fields)
    %alphaSpec = dividend_alpha ./ denominator_alpha;
    h = [];
    %h(end+1) = plot(depths,alphaSpec,'LineWidth',2.5,'LineStyle',':');
    %hold on;
    % xlabel('depths in mm')
    % ylabel('alpha in Gy^{-1}')
    % legend('LEMIII014')


    % for i = [1 2 3]
        %tableEntry.data(i) = RBEtable.(fields{i,1});
        % modelEnergiesMev_single = (RBEtable_LEMI.data(2).energies)/u; %MeV/u
        %modelEnergiesMeV_single = modelEnergies/u;
        % if i == 1
        %     idx = find(strcmp(fields, 'LEMI'));
        % elseif i == 2
        %     idx = find(strcmp(fields,'LEMII'));
        % elseif i == 3
        %     idx = find(strcmp(fields,'LEMIII'));
        % end

        alphaInterp_single = interp1(modelEnergies,RBEtable.data(6).alpha,specEnergies_single,"linear","extrap");
        betaInterp_single  = interp1(modelEnergies,RBEtable.data(6).beta,specEnergies_single,"linear","extrap");
    % end

    %for i = [1 2 3]
        % Plot first group (1,2,3)
        %if ismember(i,[3])
        figure
        dividend_alpha_table_single = dividend_alpha_zero + sum(alphaInterp_single.*dEdxInterp_rbetable_single.*spectrum(particle(6)).fluenceSpectrum',2); % sum over energy
        %dividend_alpha_table_single = dividend_alpha_zero + sum(alphaInterp_single.*spectrum(particle(6)).fluenceSpectrum',2); % sum over energy
        alphaSpec_single = (dividend_alpha_table_single ./ denominator_alpha_table_single); 
        h(end+1) = plot(depths,alphaSpec_single,'LineWidth',2.5,'LineStyle',':');
        hold on;
        h(end+1) = plot(MCFdepths,MCFalpha);
        %end

        % Plot second group (6,7)
        % if ismember(i,[6 7])
        %     alphaInterp_table(i,:) = interp1(modelEnergies,tableEntry.data(i).data(6).alpha(:,table),specEnergies,"linear","extrap");
        %     dividend_alpha_table(:,i) = dividend_alpha_zero + sum(alphaInterp_table(i,:).*dEdxInterp.*spectrum.fluenceSpectrum',2);
        %     alphaSpec(:,i) = dividend_alpha_table(:,i) ./ denominator_alpha_table;
        %     h(end+1) = plot(depths,alphaSpec(:,i),'LineWidth',2.5,'LineStyle',':');
        %     hold on;
        % end
    %end

    % Assign legend to stored handles (order matters!)
    % legend(h, {'LEMI','LEMII','Scholz2006 0.1','Scholz2006 0.14','Russo2011 0.14'});
    xlabel('depths in mm')
    ylabel('alpha in Gy^{-1}')
    legend('RBEtable MCF MKM','MCF MKM machine');

    %%

     figure
    %for i = [1 2 3]
        dividend_sqrtbeta_table_single = dividend_sqrtbeta_zero + sum(betaInterp_single.*dEdxInterp_rbetable_single.*spectrum(particle(6)).fluenceSpectrum',2);
        sqrtbetaSpec_single = (dividend_sqrtbeta_table_single ./ denominator_sqrtbeta_table_single);
        h(end+1) = plot(depths,sqrtbetaSpec_single,'LineWidth',2.5,'LineStyle',':');
        hold on;
        plot(MCFdepths,MCFbeta);
    % end
    xlabel('depths in mm')
    ylabel('beta in Gy^{-2}')
    legend('RBEtable MCF MKM','MCF MKM machine');
    % legend('LEMI','LEMII','LEMIII');


%end

% alphaSpectrum = sum(cat(3,alphaSpec{:}),3);
% betaSpectrum = sum(cat(3,sqrtbetaSpec{:}),3);
% figure
% h(end+1) = plot(depths, alphaSpec_single(:,i),'LineWidth',2.5,'LineStyle',':');
% hold on
%     xlabel('depths in mm')
%     ylabel('alpha in Gy^{-1}')
%     legend('LEMI','LEMII','LEMIII');
% figure
% h(end+1) = plot(depths,betaSpectrum(:,3),'LineWidth',2.5,'LineStyle',':');
% hold on
%     xlabel('depths in mm')
%     ylabel('beta in Gy^{-2}')
%     legend('LEMI','LEMII','LEMIII');


%% Sum up the secondaries and primaries
%figure
% fields = fieldnames(RBEtable);

MCFdepths = MCF_machine.data(76).depths;
MCFalpha = MCF_machine.data(76).alpha;
MCFbeta = MCF_machine.data(76).beta;

strings = ["H","He","Li","Be","B","C"];
As = [1 4 7 9 10 12];
particle = [1 4 5 6 7 8];

bdEntry = machine.data(48);


depths = bdEntry.depths;
depthGrad = gradient(depths);
fluence = machine.data(48).Fluence;

energies = fluence.spectra(8).energyBin;
spectrum = fluence.spectra;


dividend_alpha_zero = zeros(size(depths));
denominator_alpha_zero = zeros(size(depths));
dividend_sqrtbeta_zero = zeros(size(depths));
denominator_sqrtbeta_zero = zeros(size(depths));

specEnergies = cell(1,6);
dEdxEnergy = cell(1,6);
dEdx = cell(1,6);
dEdxInterp_rbetable = cell(1,6);
alphaInterp = cell(1,6);
betaInterp = cell(1,6);
dividend_alpha_table = cell(1,6);
dividend_sqrtbeta_table = cell(1,6);
alphaSpec = cell(1,6);
sqrtbetaSpec = cell(1,6);


for p = 1:6
    modelEnergies = RBEtable.data(p).energies; % MeV/u
    currZ = spectrum(particle(p)).Z;

    %modelEnergiesMev(:,p) = modelEnergies/u;
    dEdxEnergy{p} = SPtable.data(p).energies;
    dEdx{p} = SPtable.data(p).dEdx;
    specEnergies{p} = spectrum(particle(p)).energyBin; %MeV/u
    dEdxInterp_rbetable{p} = interp1(dEdxEnergy{p},dEdx{p}(:),specEnergies{p},"linear","extrap");
    denominator_alpha_table(:,p) = denominator_alpha_zero + sum(dEdxInterp_rbetable{p}.*spectrum(particle(p)).fluenceSpectrum',2); %sum over energy
    denominator_sqrtbeta_table(:,p) = denominator_sqrtbeta_zero + sum(dEdxInterp_rbetable{p}.*spectrum(particle(p)).fluenceSpectrum',2);

    %for i = 1:numel(fields)
    %alphaSpec = dividend_alpha ./ denominator_alpha;
    h = [];
    %h(end+1) = plot(depths,alphaSpec,'LineWidth',2.5,'LineStyle',':');
    %hold on;
    % xlabel('depths in mm')
    % ylabel('alpha in Gy^{-1}')
    % legend('LEMIII014')


    % for i = [1 2 3]
    %     tableEntry.data(i) = RBEtable.(fields{i,1});
    %     modelEnergiesMev(:,p) = (RBEtable.(fields{i,1}).data(p).energies)/u;
    %     if i == 1
    %         idx = find(strcmp(fields, 'LEMI'));
    %     elseif i == 2
    %         idx = find(strcmp(fields,'LEMII'));
    %     elseif i == 3
    %         idx = find(strcmp(fields,'LEMIII'));
    %     end
        alphaInterp{p} = interp1(modelEnergies,RBEtable.data(p).alpha(:,1),specEnergies{p},"linear","extrap");
        betaInterp{p}  = interp1(modelEnergies,RBEtable.data(p).beta(:,1),specEnergies{p},"linear","extrap");
    % end

    % for i = [1 2 3]
        % Plot first group (1,2,3)
        %if ismember(i,[3])

        dividend_alpha_table{p} = dividend_alpha_zero + sum(alphaInterp{p}.*dEdxInterp_rbetable{p}.*spectrum(particle(p)).fluenceSpectrum',2); % sum over energy
        %alphaSpec{p} = (dividend_alpha_table{p} ./ denominator_alpha_table(:,p));
        % h(end+1) = plot(depths,alphaSpec{p}(:,i),'LineWidth',2.5,'LineStyle',':');
        % hold on;
        %end

        % Plot second group (6,7)
        % if ismember(i,[6 7])
        %     alphaInterp_table(i,:) = interp1(modelEnergies,tableEntry.data(i).data(6).alpha(:,table),specEnergies,"linear","extrap");
        %     dividend_alpha_table(:,i) = dividend_alpha_zero + sum(alphaInterp_table(i,:).*dEdxInterp.*spectrum.fluenceSpectrum',2);
        %     alphaSpec(:,i) = dividend_alpha_table(:,i) ./ denominator_alpha_table;
        %     h(end+1) = plot(depths,alphaSpec(:,i),'LineWidth',2.5,'LineStyle',':');
        %     hold on;
        % end
    %end

    % Assign legend to stored handles (order matters!)
    % legend(h, {'LEMI','LEMII','Scholz2006 0.1','Scholz2006 0.14','Russo2011 0.14'});
    % xlabel('depths in mm')
    % ylabel('alpha in Gy^{-1}')
    % legend('LEMI','LEMII','LEMIII');

    % figure
    % for i = [1 2 3]
        dividend_sqrtbeta_table{p} = dividend_sqrtbeta_zero + sum(betaInterp{p}.*dEdxInterp_rbetable{p}.*spectrum(particle(p)).fluenceSpectrum',2);
        %sqrtbetaSpec{p} = (dividend_sqrtbeta_table{p} ./ denominator_sqrtbeta_table(:,p));
        % h(end+1) = plot(depths,sqrtbetaSpec(:,i),'LineWidth',2.5,'LineStyle',':');
        % hold on;
    % end
    % 
    % xlabel('depths in mm')
    % ylabel('beta in Gy^{-2}')
    % legend('LEMI','LEMII','LEMIII')

end

dividend_alpha_sum = sum(cat(3,dividend_alpha_table{:}),3);
denominator_alpha_sum = sum(denominator_alpha_table,2);
total_alphaSpec = dividend_alpha_sum ./ denominator_alpha_sum;
dividend_beta_sum = sum(cat(3,dividend_sqrtbeta_table{:}),3);
denominator_beta_sum = sum(denominator_sqrtbeta_table,2);
total_betaSpec = dividend_beta_sum ./ denominator_beta_sum;

figure
h(end+1) = plot(depths,total_alphaSpec,'LineWidth',2.5,'LineStyle',':');
hold on 
h(end+1) = plot(MCFdepths,MCFalpha);

figure
h(end+1) = plot(depths,total_betaSpec,'LineWidth',2.5,'LineStyle',':');
hold on 
h(end+1) = plot(MCFdepths,MCFbeta);

% betaSpectrum = sum(cat(3,sqrtbetaSpec{:}),3);
%%
figure
h(end+1) = plot(depths,alphaSpectrum(:,1),'LineWidth',2.5,'LineStyle',':');
hold on
    xlabel('depths in mm')
    ylabel('alpha in Gy^{-1}')
    legend('LEMI','LEMII','LEMIII');
figure
h(end+1) = plot(depths,betaSpectrum(:,3),'LineWidth',2.5,'LineStyle',':');
hold on
    xlabel('depths in mm')
    ylabel('beta in Gy^{-2}')
    legend('LEMI','LEMII','LEMIII');

% alpha0 = [0.0708, 0.313, 0.228]; %0.282;
% beta0 = [0.0615, 0.0615, 0.02];
% % alphaX = tableEntry.alphaX;
% % betaX  = tableEntry.betaX;
% 
% for j = 1:size(alpha0,2)
%     alphaMKM(:,j) = (zMix*beta0(j) + alpha0(j));
%     plot(depths,alphaMKM(:,j),'LineWidth',2.5);
% end

%% Interpolation

sizeAlphaSingle = (1:length(alphaSpec_single))';
sizeNEWSpec = linspace(1,length(alphaSpec_single),length(machine.data(76).alpha))';

alphaInterp_MCF = interp1(sizeAlphaSingle,alphaSpec_single,sizeNEWSpec,'spline');
% 'spline' oder 'pchip' für glattere Kurven

SizeAlphaComp = size(alphaInterp_MCF);

Ratio = alphaInterp_MCF ./ machine.data(76).alpha;

x_ratio = (1:length(Ratio))';
x_target = linspace(1,length(Ratio),length(alphaSpec_single))';

ratio_interp = interp1(x_ratio,Ratio,x_target,'spline');

result = alphaSpec_single ./ ratio_interp;
result2 = alphaSpec_single ./  1.9125;

x_RBEtable = (1:length(RBEtable.data(6).alpha))';
x_RBEtarget = linspace(1,length(RBEtable.data(6).alpha),length(alphaSpec_single))';

RBEtableAlpha_interp = interp1(x_RBEtable,RBEtable.data(6).alpha,x_RBEtarget,'spline');


result3 = alphaSpec_single ./ RBEtableAlpha_interp;

