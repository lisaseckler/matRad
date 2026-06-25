f = figure('WindowState','maximized');

%HIT = machine;

blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
purple = [0.4940 0.1840 0.5560];
brown = [0.6350 0.0780 0.1840];
colors = [blue; blue; blue; orange; green; red; purple; brown];

lineWidth = 2;
markerSize = 15;

x = HIT.data(95).depths;

vois = [1 4 5 6 7 8];

for i = vois
    fluenceHIT = sum(HIT.data(95).Fluence.spectra(i).fluenceSpectrum,1);
    plot(x,fluenceHIT,'-', 'DisplayName',num2str(i),'Color', colors(i,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
    c =1;
    leg = {};
    names = {};
    %custom legend
    for i = vois
        leg{c} = plot(nan,'Color',colors(i,:),'LineWidth',2);
        hold on
        names{c} = num2str(i);
        c = c+1;

    end
end

%MCF = machine;
x = MCF.data(72).depths;

for i = [1 4 5 6 7 8]
    fluenceMCF = sum(MCF.data(72).Fluence.spectra(i).fluenceSpectrum,1);
    plot(x,fluenceMCF,'--', 'DisplayName',num2str(i),'Color', colors(i,:), 'LineWidth', lineWidth, 'MarkerSize',markerSize)
    hold on
end


grid on;
xlabel('Depth [mm]', 'FontSize',20);
ylabel('Fluence spectra','FontSize',20);
yscale("log")
axis([0,180,0,100])
ax = gca();

leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-','LineWidth',2);
names{c} = 'HIT';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','--','LineWidth',2);
names{c} = 'MCF';
c = c+1;

legend ([leg{:}],names ) 


%% Normalisation

NormFactor = [];
%newMCF = machine;

for i = 1:166
    sumFluence = sum(MCF.data(i).Fluence.spectra(8).fluenceSpectrum,1);
    maximum = max(sumFluence(:));
    NormFactor = [NormFactor maximum];
end
%%
for i = 1:166
    for j = 1:10
        MCF.data(i).Fluence.spectra(j).fluenceSpectrum = MCF.data(i).Fluence.spectra(j).fluenceSpectrum ./ NormFactor(j);
    end
end

%% saving
MCF.meta.machine = 'MCF_latestVersion_correctedFluence_Americans2';
machine = MCF;

save("carbon_MCF_latestVersion_correctedFluence_Americans2.mat","machine")

%% setting the first energyBins from MCF to the ones from HIT

%energyBins = [HIT.data(14).Fluence.spectra(1).energyBin(1) HIT.data(14).Fluence.spectra(2).energyBin(1) HIT.data(14).Fluence.spectra(3).energyBin(1) HIT.data(14).Fluence.spectra(4).energyBin(1) HIT.data(14).Fluence.spectra(5).energyBin(1) HIT.data(14).Fluence.spectra(6).energyBin(1) HIT.data(14).Fluence.spectra(7).energyBin(1) HIT.data(14).Fluence.spectra(8).energyBin(1) HIT.data(14).Fluence.spectra(9).energyBin(1) HIT.data(14).Fluence.spectra(10).energyBin(1) HIT.data(14).Fluence.spectra(11).energyBin(1)];

for i = 1:166
    for j = 5:7
        % n = length(MCF.data(i).Fluence.spectra(j).fluenceSpectrum(:,2));
        % v1 = HIT.data(1).Fluence.spectra(1).fluenceSpectrum(:,2);  % Spaltenvektor
        % 
        % result = [v1; zeros(n - length(v1), 1)];  % vertikale Konkatenation mit ;
        % result = result(1:n);
        % 
        % MCF.data(i).Fluence.spectra(j).fluenceSpectrum(:,2) = result;
        machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,1) = 0;
        machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,2) = 0;
        machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,3) = 0;
        machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,4) = 0;
        machine.data(i).Fluence.spectra(j).fluenceSpectrum(:,5) = 0;
    end
end

