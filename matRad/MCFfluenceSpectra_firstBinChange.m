%% Changing MCF fluence spectra for the first energyBin
clear
% loading machines
HIT = load('basedata\carbon_HITfluenceSpectra_updatedBig.mat');
MCF = load('basedata\carbon_MCF_Fluence_latestVersion.mat');

%% calculating correction factor
fraIx = 1;
fragmentSum = {};
depthSum = [];

for i = 1:166
    for j = [1 4 5 6 7]
        for k = 1: numel(MCF.machine.data(i).depths)
            initSum = trapz(MCF.machine.data(i).Fluence.spectra(j).energyBin,MCF.machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,k));
            depthSum = [depthSum initSum];
        end
        fragmentSum{i,fraIx} = depthSum;
        fraIx = fraIx + 1;
        depthSum = [];
    end
    fraIx = 1;
end


%% setting the first bin to the second

for i = 1:166
    for j = [1 4 5 6 7]
        MCF.machine.data(i).Fluence.spectra(j).fluenceSpectrum(1,:) = MCF.machine.data(i).Fluence.spectra(j).fluenceSpectrum(2,:);
    end
end

% calculating the integral sum again
corrIx = 1;
corrSum = {};
depthSum = [];

for i = 1:166
    for j = [1 4 5 6 7]
        for k = 1: numel(MCF.machine.data(i).depths)
            initSum = trapz(MCF.machine.data(i).Fluence.spectra(j).energyBin,MCF.machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,k));
            depthSum = [depthSum initSum];
        end
        corrSum{i,corrIx} = depthSum;
        corrIx = corrIx + 1;
        depthSum = [];
    end
    corrIx = 1;
end

%% correcting the first bin
corrTerm = {};

for i = 1:166
    for j = 1:5
        corrTerm{i,j} = corrSum{i,j} ./ fragmentSum{i,j};
    end
end

%% correcting the first bin in the fluence spectra

fragmentIx = 1;

for i = 1:166
    for j = [1 4 5 6 7]
        for k = 1: numel(MCF.machine.data(i).depths)
            MCF.machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,k) = MCF.machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,k) ./ corrTerm{i,fragmentIx}(1,k);
        end
        fragmentIx = fragmentIx + 1;
    end
    fragmentIx = 1;
end

%% plotting to compare HIT and MCF

figure;
plot(HIT.machine.data(95).Fluence.spectra(4).energyBin,HIT.machine.data(95).Fluence.spectra(4).fluenceSpectrum(:,8));
hold on;
plot(MCF.machine.data(72).Fluence.spectra(4).energyBin,MCF.machine.data(72).Fluence.spectra(4).fluenceSpectrum(:,6));

xlabel('energyBin')
ylabel('fluenceSpectrum helium at 112.5 mm peak position in 13 mm depth')

%% saving new MCF machine

machine = MCF.machine;

machine.meta.machine = 'MCF_latestVersion_correctedFirstBin';

save("carbon_MCF_latestVersion_correctedFirstBin","machine")