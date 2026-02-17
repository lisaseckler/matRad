%% extract MCF-MKM table from excel

filename = 'Alpha_Beta_Tables_MCFMKM.xlsx';

opts = detectImportOptions(filename);
opts.VariableNamesRange = 2;   % headers in row 2
opts.DataRange = 3;            % data starts in row 3

T = readtable(filename, opts);

%% Check actual variable names
T.Properties.VariableNames

%% Find unique (Z, A) pairs
[ZA_unique, ~, idx] = unique(T(:, {'Z','A'}), 'rows');

nZA = height(ZA_unique);

% Preallocate output table
Out = table;
Out.Z        = ZA_unique.Z;
Out.A        = ZA_unique.A;
Out.energies = cell(nZA,1);
Out.alpha    = cell(nZA,1);
Out.beta     = cell(nZA,1);

% Collect all energies/alpha/beta for each (Z, A)
for i = 1:nZA
    rows = (idx == i);

    Out.energies{i} = T.KE_MeV_u_(rows);
    Out.alpha{i}    = T.Alpha_Gy__1_(rows);
    Out.beta{i}     = T.Beta_Gy__2_(rows);
end

%% restructure
MCF = Out;
load("RBEtable_LEMI_Scholz06_AX313_BX615.mat")

MCFMKM = RBEtable;
%%

MCFMKM.meta.model = 'MCFMKM';
MCFMKM.meta.description = 'Model created from Shannon Hartzell';
MCFMKM.meta.modelParameters.alphaX = 0.117;
MCFMKM.meta.modelParameters.betaX = 0.0615;
MCFMKM.meta.modelParameters.rNucleus = 4.5;
MCFMKM.meta.modelParameters.rDomain = 0.28;

% delete calculusType

for i = 1:8
    MCFMKM.data(i).energies = MCF.data.energies{i};
    MCFMKM.data(i).alpha = MCF.data.alpha{i};
    MCFMKM.data(i).beta = MCF.data.beta{i};
end

RBEtable = MCFMKM;

%%
RBEtable.data(7).Z = 7;
RBEtable.data(8).Z = 8;
RBEtable.data(7).A = 14;
RBEtable.data(8).A = 16;

%%
save("RBEtable_MCFMKM_photon_Shannon_AX117_BX615","RBEtable")