function fitData = fitCarbonBasedata(EnergyDepositedZR, Data_Fluence, params, zdepth, rdepth, scoringArea, visBool)


%fileID = fopen(fullfile('E:\Code\CreateHeliumBasedata\SimulationHelium\fitHeliumBasedata.log'),'a');
%fprintf(fileID, 'Starting Simu %g \n', params.energy);

nZspec = numel([Data_Fluence.spectra(:).Z]);
N = size(EnergyDepositedZR,1);

% extract IDD 
IDD = sum(EnergyDepositedZR,2);
% interpolate range at 80% dose after peak.
% Because of Helium High Dose Max Energy  change fitting start value for
% max
range_estimated = params.estRange;
[~, r_min] = min(abs(zdepth - (range_estimated - 30)));
[~, r_max] = min(abs(zdepth - (range_estimated + 30)));
[maxIDD, max_ix] = max(IDD(r_min:r_max));
max_ix = r_min + max_ix -1;
[~, r80ix] = min(abs(IDD(max_ix:end) - 0.8 * maxIDD));
r80ix = r80ix-1; %max_ix index is 1
r80 = interp1(IDD(max_ix + r80ix - 1:max_ix + r80ix + 1), ...
                 zdepth(max_ix + r80ix - 1:max_ix + r80ix + 1), 0.8 * maxIDD);

% fit gaussians to radial Dose slices to get sigma 
sigma = zeros(N,2);
weight = zeros(N,1);

% fit sigmas till dose reaches, a percent of max dose
[~, ixSigmaFit] = min(abs(IDD(max_ix:end) - maxIDD));
if zdepth(max_ix+ixSigmaFit-1) < r80 + 20
    [~, ixSigmaFit] = min(abs(zdepth(max_ix:end) - (r80+20)));
end

% Fitting Functions
gauss1       = @(p,x) 1./(2*pi*p.^2) * exp(-x.^2./(2*p.^2));
gauss2       = @(w1, p1, p2, x) (1-w1) * gauss1(p1,x) + w1 * gauss1(p2,x);
gauss3       = @(w1, w2, p1, p2, p3, x) (1-w1-w2) * gauss1(p1,x) + w1 * gauss1(p2,x) + w2*gauss1(p3,x);
objFunc    = fittype(gauss2);
objFuncSpec = fittype(gauss3);


for i = 1:N %max_ix+ixSigmaFit-1

    % normaliseProfile
    profile = EnergyDepositedZR(i,:);
    profile = profile / (sum(profile)); 
    profile = profile./scoringArea;
   
    
   % double Gaußfit
    start = [0.1, params.sigma0+1,  params.sigma0+5];
    lb = [0.001, params.sigma0, params.sigma0];
    ub = [1, params.sigma0+100 params.sigma0+100];

    weightFit = zeros(size(profile'));
    weightFit(profile>0) = 1./profile(profile>0);


    [fitObject,gof,output] = fit(rdepth', profile', objFunc, ...
        'Start', [start], 'Lower', [lb], 'Upper', [ub],  'Weight', weightFit);

    if output.exitflag < 1
        msg = sprintf('Fitting of double gauß for energy %g and slice %g at depth %g mm stopped with exitflag %g', params.energy, i, zdepth(i), output.exitflag);
        msg = [msg,' Error Message : ', output2.message, '\n'];
        msg = convertCharsToStrings(msg);
        warning(msg);
    end

    [sigma(i,:), ixSortedSigma] = sort([fitObject.p1, fitObject.p2]);
    wTmp = [1-fitObject.w1,fitObject.w1];
    wTmp = wTmp(ixSortedSigma(2:end));
    if sum(wTmp) >= 1
        warning('Wrong weights fited for energy %g and slice %g \n', params.energy, i);
    end
    weight(i) = wTmp;

    if visBool && ((i==1) || (mod(i,25) == 0))
        figure()
        sgtitle(convertCharsToStrings(sprintf('Energy %g MeV/u at depths %g mm', params.energy, zdepth(i))))
       
        
        subplot(1,2,1)
        plot(rdepth, profile, '.')
        hold on
        plot(rdepth, gauss2(fitObject.w1,fitObject.p1, fitObject.p2, rdepth))
        xlabel('Lateral Depth [mm]')
        ylabel('Profile [MeV /mm^3]')
        title('Gauss2 Fit')

        subplot(1,2,2)
        semilogy(rdepth, profile, '.')
        hold on
        semilogy(rdepth, gauss2(fitObject.w1,fitObject.p1, fitObject.p2, rdepth))
        xlabel('Lateral Depth [mm]')
        ylabel('Profile [MeV /mm^3]')
        title('Gauss2 Fit')
        
        fprintf('Done with Energy %g MeV/u at depths %g mm \n', params.energy, zdepth(i))
    end
    
    
end



% only inlcude where sigma whas fited and correct for inital sigma, spot size
zdepth  = zdepth; %(1:max_ix +ixSigmaFit-1);
IDD     = IDD; %(1:max_ix +ixSigmaFit-1);
weight = weight; %(1:max_ix +ixSigmaFit-1,:);
%sigma  = sqrt(sigma(1:max_ix +ixSigmaFit-1,:).^2 - params.sigma0.^2);
sigma = sqrt(sigma.^2 - params.sigma0.^2);

% for ixZ = 1:nZspec
%     Data_Fluence.spectra(ixZ).fluenceSpectrum = Data_Fluence.spectra(ixZ).fluenceSpectrum(:,1:max_ix +ixSigmaFit-1);
% end


if ~isreal(sigma)
    msg = sprintf('sigma of energy %g contains complex values, CHECK !!!', params.energy);
    msg = convertCharsToStrings(msg);
    warning(msg);
end

% conversion factor for IDDs: convert to MeV * cm^2/g per primary
% input MeV/mm
% Divide Mass of Water rho = 1 g/cm^3
cf = 10;

% Add 0 depth
% interpolate sigma on depths of IDD
IDD     = interp1(zdepth, IDD, [0,zdepth], 'linear', 'extrap')';
weight = interp1(zdepth, weight, [0,zdepth], 'linear', 'extrap');
sigma   = [[0,0];sigma];

[XE,YZ] = meshgrid(Data_Fluence.energyBin,[0,zdepth]);
[xe,yz] = meshgrid(Data_Fluence.energyBin,zdepth);
ANorm = [1,2,3,4,7,9,11,12,14,16];
for ixZ = 1:nZspec-1
    Data_Fluence.spectra(ixZ).fluenceSpectrum = interp2(xe,yz, Data_Fluence.spectra(ixZ).fluenceSpectrum',XE,YZ, 'linear')';
    Data_Fluence.spectra(ixZ).energyBin = Data_Fluence.energyBin./ANorm(ixZ);
end
ixZ = ixZ+1;
[XE,YZ] = meshgrid(Data_Fluence.energyBinEl,[0,zdepth]);
[xe,yz] = meshgrid(Data_Fluence.energyBinEl,zdepth);
Data_Fluence.spectra(ixZ).fluenceSpectrum = interp2(xe,yz, Data_Fluence.spectra(ixZ).fluenceSpectrum',XE,YZ, 'linear')';
Data_Fluence.spectra(ixZ).energyBin = Data_Fluence.energyBinEl;

for ixZ = 1:nZspec
    Data_Fluence.spectra(ixZ).fluenceDepth = sum(Data_Fluence.spectra(ixZ).fluenceSpectrum,1);
end

zdepth  = [0,zdepth];
%% get rid of possible negative values thru interpolation
IDD(IDD<0)          = 0;
sigma(sigma<0)      = 0;
weight(weight<0)    = 0;

% save data in machine
fitData.energy      = params.energy;
fitData.range       = r80;    
fitData.depths      = zdepth';
fitData.Z           = IDD* cf;
fitData.peakPos     = zdepth(max_ix+1); %zeroDepth added
fitData.sigma       = sigma;
fitData.weight      = weight;
fitData.Fluence.energyBin = Data_Fluence.energyBin;
fitData.Fluence.energyBinEl = Data_Fluence.energyBinEl;
fitData.Fluence.spectra = Data_Fluence.spectra;

%fprintf(fileID, 'Finished Simu %g \n', params.energy);
%fclose(fileID);

end

