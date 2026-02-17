%% Interpolation sigma MCF and Simonas base data

% Energies from the base data
EnergiesSimona = arrayfun(@(x) x.energy, machine.data);  % 121x1
EnergiesMCF = arrayfun(@(x) x.energy, MCF.data);  % 166x1

numberEnergiesMCF = numel(MCF.data);
idxLow  = zeros(numberEnergiesMCF,1);
idxHigh = zeros(numberEnergiesMCF,1);
w       = zeros(numberEnergiesMCF,1);

for j = 1:numberEnergiesMCF
    if EnergiesMCF(j) <= EnergiesSimona(1)
        idxLow(j) = 1;
        idxHigh(j) = 1;
        w(j) = 0;
    elseif EnergiesMCF(j) >= EnergiesSimona(end)
        idxLow(j) = numel(EnergiesSimona);
        idxHigh(j) = numel(EnergiesSimona);
        w(j) = 0;
    else
        idxHigh(j) = find(EnergiesSimona >= EnergiesMCF(j), 1, 'first');
        idxLow(j)  = idxHigh(j) - 1;
        w(j) = (EnergiesMCF(j) - EnergiesSimona(idxLow(j))) / (EnergiesSimona(idxHigh(j)) - EnergiesSimona(idxLow(j)));
    end
end

diff = idxHigh - idxLow;
%%
for i = 1:numberEnergiesMCF
    if diff(i) == 0
        otherwayAround(i) = 1;
    else
        otherwayAround(i) = 0;
    end
    
end

% here are the energys that needs to be extrapolated
energyIdx = find(otherwayAround);
%%
for k = 1:numel(EnergiesSimona)
    if EnergiesSimona(k) <= EnergiesMCF(k+8)
        diffEnergies = abs(EnergiesMCF(k+8) - EnergiesSimona(k));
        diffEnergies2 = abs(EnergiesMCF(k+8) - EnergiesSimona(k+1));
        if diffEnergies > diffEnergies2
            idx(k) = k+1;
        elseif diffEnergies2 > diffEnergies
            idx(k) = k;
        end
    elseif EnergiesSimona(k) >= EnergiesMCF(k+8)
        diffEnergies = abs(EnergiesSimona(k) - EnergiesMCF(k+8));
        diffEnergies2 = abs(EnergiesMCF(k+8) - EnergiesSimona(k-1));
        if diffEnergies > diffEnergies2
            idx(k) = k-1;
        elseif diffEnergies2 > diffEnergies
            idx(k) = k;
        end
    end
end

%% try with while loop

l = 1;
for k = 1:numel(EnergiesMCF)
    if EnergiesSimona(l) <= EnergiesMCF(k)
        if l <= 120
            diffEnergies = abs(EnergiesMCF(k) - EnergiesSimona(l));
            diffEnergies2 = abs(EnergiesMCF(k) - EnergiesSimona(l+1));
            if diffEnergies > diffEnergies2
                idx(k) = l+1;
                l = l+1;
            elseif diffEnergies2 > diffEnergies
                idx(k) = l;
            end
        elseif l == 121
            idx(k) = l;
        end
    elseif EnergiesSimona(l) >= EnergiesMCF(k)
        if l > 1
            diffEnergies = abs(EnergiesSimona(l) - EnergiesMCF(k));
            diffEnergies2 = abs(EnergiesMCF(k) - EnergiesSimona(l-1));
            if diffEnergies > diffEnergies2
                idx(k) = l-1;
            elseif diffEnergies2 > diffEnergies
                idx(k) = l;
            end
        elseif l == 1
            idx(k) = l;
        end
    end
end

%%
% yq = interp1(x,y,xq)

for i = 1:166
    MCF.data(i).sigma = interp1(machine.data(idx(i)).depths/machine.data(idx(i)).peakPos,machine.data(idx(i)).sigma,MCF.data(i).depths/MCF.data(i).peakPos);
    MCF.data(i).sigma1 = interp1(machine.data(idx(i)).depths/machine.data(idx(i)).peakPos,machine.data(idx(i)).sigma1,MCF.data(i).depths/MCF.data(i).peakPos);
    MCF.data(i).sigma2 = interp1(machine.data(idx(i)).depths/machine.data(idx(i)).peakPos,machine.data(idx(i)).sigma2,MCF.data(i).depths/MCF.data(i).peakPos);
    MCF.data(i).sigma1 = interp1(machine.data(idx(i)).depths/machine.data(idx(i)).peakPos,machine.data(idx(i)).sigma1,MCF.data(i).depths/MCF.data(i).peakPos);
    MCF.data(i).weight = interp1(machine.data(idx(i)).depths/machine.data(idx(i)).peakPos,machine.data(idx(i)).weight,MCF.data(i).depths/MCF.data(i).peakPos);
end

for i = 1:166
    MCF.data(i).sigma(1) = 0;
    MCF.data(i).sigma1(1) = 0;
end

MCF.meta.dataType = 'doubleGauss';

%% saving
machine = MCF;
save("carbon_MCF_Fluence_updated_withoutRBE.mat")

load("carbon_MCF_Fluence_updated.mat")

machine.meta.dataType = 'doubleGauss';

for i = 1:166
    machine.data(i).sigma = MCF.data(i).sigma;
    machine.data(i).sigma1 = MCF.data(i).sigma1;
    machine.data(i).sigma2 = MCF.data(i).sigma2;
    machine.data(i).weight = MCF.data(i).weight;    
end

save("carbon_MCF_Fluence_updated.mat","machine")