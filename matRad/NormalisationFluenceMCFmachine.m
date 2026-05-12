f = figure('WindowState','maximized');

%machine_HIT = machine;

blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
purple = [0.4940 0.1840 0.5560];
brown = [0.6350 0.0780 0.1840];
colors = [blue; blue; blue; orange; green; red; purple; brown];

lineWidth = 2;
markerSize = 15;

x = machine_HIT.data(95).depths;

vois = [1 4 5 6 7 8];

for i = vois
    fluenceHIT = sum(machine_HIT.data(95).Fluence.spectra(i).fluenceSpectrum,1);
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

%newMCF = machine;
x = newMCF.data(72).depths;

for i = [1 4 5 6 7 8]
    fluenceMCF = sum(newMCF.data(72).Fluence(i).fluenceSpectrum,1);
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
newMCF = machine;

for i = 1:166
    sumFluence = sum(newMCF.data(i).Fluence(8).fluenceSpectrum,1);
    maximum = max(sumFluence(:));
    NormFactor = [NormFactor maximum];
end
%%
for i = 1:166
    for j = 1:8
        newMCF.data(i).Fluence(j).fluenceSpectrum = newMCF.data(i).Fluence(j).fluenceSpectrum ./ NormFactor(j);
    end
end

%% saving
machine_MCF.meta.machine = 'MCF_Fluence_updated_withoutRBE_origSingleGauss_updatedNormFluence';
machine = machine_MCF;

save("carbon_MCF_Fluence_updated_withoutRBE_origSingleGauss_updatedNormFluence.mat","machine")