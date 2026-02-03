%% Introduction
% This script loads the TOPAS simulation results for all energys, performes
% the fit for each energy and writes it all into one nice machine basedata
% file
% Unit of Length mm
% Unit of Energy MeV/u if not staede otherwise

 load("SimuParams.mat")
 load("BeamOptics.mat")
 machine_Carbon = load("carbon_HITfixedBL.mat");


%% Init Machine

machine = struct;

machine.meta = struct;
machine.meta.radiationMode = 'carbon';
machine.meta.dataType = 'multipleGauss';
machine.meta.created_on = date;
machine.meta.created_by = 'matRad Team';
machine.meta.description = 'carbon beamdata for HIT machine, simulated with TOPAS';
machine.meta.SAD = machine_Carbon.machine.meta.SAD;
machine.meta.BAMStoIsoDist = machine_Carbon.machine.meta.BAMStoIsoDist;
machine.meta.machine = 'HITfixedBL_TOPAS';
machine.meta.LUTspotSize = machine_Carbon.machine.meta.LUT_bxWidthminFWHM;

machine.data = struct;

params = struct();

params.SAD = machine.meta.SAD;
params.r = 1;
params.rNum = 70;
params.zCoarse = 2;
params.zFine = 0.2;
params.zFineNum = 100;
params.numPrimary = 1000000;
params.Z = 6;
params.A = 12;

dataType = 'binary';

%% Fitting Loop

foldername = 'C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\matRad\Scripts_big\';

visBool = false;
for i = 113

    if ismember(i, [113])
        visBool = true;
    else
        visBool = false;
    end

    % Parameters for Fitting;
    params.estRange = approxRange(i);
    params.sigma0 = sigmaSqIso(i,:);
    params.energy = E(i);
    params.zCoarseNum = simulateHalfLength(i);
    params.spectrum_ebinsNum = Ebins(i);
    params.spectrum_ebinsMax = Emax(i);
    params.electron_spectrum_ebinsNum = 100;
    params.electron_spectrum_ebinsMax=1;
    % Beam energy
    params.BeamOptics.BeamEnergyBins = EnergyBins{i};
    params.BeamOptics.BeamEnergySpectrum = EnergySpectrum{i};
    params.BeamOptics.SigmaX = sigmaSqIso(i,:);
    params.BeamOptics.SigmaY = params.BeamOptics.SigmaX;
    params.BeamOptics.SigmaXprime = sigmaT(i,:);
    params.BeamOptics.SigmaYprime = params.BeamOptics.SigmaXprime;
    params.BeamOptics.CorrelationX = rho(i,:);
    params.BeamOptics.CorrelationY = params.BeamOptics.CorrelationX;

    %  Data Dose and RBE
    Data_Dose =  loadSimuDataDose(foldername, params, E(i),'binary' ); 
    Data_Fluence = loadSimuDataFluence(foldername, params, E(i),'binary' ); 

    % Dont normalize to rBin width because for lataral fits i will
    %normalize to the area of the scoring rings anyway
    % Combine Scoring Volumes
    [Data_Dose, Data_Fluence, Zdepths, rlengths, scoringArea]  = combineScoringVolumes(Data_Dose, Data_Fluence, params, approxRange(i));

    % TODO fluence fit + readin
    fitData = fitCarbonBasedata(Data_Dose.ED, Data_Fluence,params, Zdepths, rlengths, scoringArea, visBool);
    machine.data(i).range                   = fitData.range;
    machine.data(i).energy                  = fitData.energy;
    machine.data(i).depths                  = fitData.depths;
    machine.data(i).Z                       = fitData.Z;
    machine.data(i).peakPos                 = fitData.peakPos;
    machine.data(i).sigmaMulti                   = fitData.sigma;
    machine.data(i).weightMulti                  = fitData.weight';
    machine.data(i).offset                  = 0; 
    machine.data(i).BeamOpticsSimulated     = params.BeamOptics;
    machine.data(i).Fluence                 = fitData.Fluence;


    % TODO
    machine.data(i).initFocus.dist = machine_Carbon.machine.data(i).initFocus.dist;
    machine.data(i).initFocus.sigma = machine_Carbon.machine.data(i).initFocus.sigma;
    machine.data(i).initFocus.SisFWHMAtIso = 2.355 .* sigmaSqIso(i,:);
    emmitance.type = 'bigaussian';
    emmitance.sigmaX    = sigmaSqIso(i,:);
    emmitance.divX      = sigmaT(i,:);
    emmitance.corrX     = rho(i,:);
    emmitance.sigmaY    = emmitance.sigmaX;
    emmitance.divY      = emmitance.divX;
    emmitance.corrY     = emmitance.corrX;
    machine.data(i).initFocus.emittance = emmitance;

    fprintf('Energy %g done \n', params.energy)
    
end


%% Save machine Data
save('carbon_HITfixedBL_TOPAS.mat', 'machine', '-v7.3')
disp('Saving Done')
