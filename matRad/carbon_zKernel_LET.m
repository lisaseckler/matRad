%% LETd comparison
clear
%% 

Zs = [1 2 3 4 5 6];
strings = ["H","He","Li","Be","B","C"]; %
As = [1 4 7 9 11 12]; %

for j = Zs
    tableEntry = RBEtable.data(Zs);
    dEdx = SPtable.data(Zs);
end

bdEntry = HIT.data(61); % HIT machine peak position 78.5 mm (61); MCF (50)


depths_HIT = bdEntry.depths;
depthGrad = gradient(depths_HIT);
fluence = bdEntry.Fluence;

fragments = [1 4 5 6 7 8]; %1 4 5 6 7 
fluenceSpectrum = fluence.spectra(fragments);

% for j = fragments
%     energies(j,:) = fluence.spectra(j).energyBin;
% end

%% plotting
figure;
hDeDx = subplot(2,2,1);
hZstar = subplot(2,2,2);
dividend_z_HIT = zeros(size(depths_HIT));
denominator_z_HIT = zeros(size(depths_HIT));
dividend_alpha = zeros(size(depths_HIT));
denominator_alpha = zeros(size(depths_HIT));
dividend_sqrtbeta = zeros(size(depths_HIT));
denominator_sqrtbeta = zeros(size(depths_HIT));
dividend_LET_HIT = zeros(size(depths_HIT));
denominator_LET_HIT = zeros(size(depths_HIT));

for i = 1:numel(fluenceSpectrum)
    %get z star
    modelEnergies = tableEntry(i).energies;
    spectrum = fluenceSpectrum(i);
    currZ = spectrum.Z;
    if currZ > 0 && currZ <= numel(As)
        rowIndices = matches(strings(Ztable.data(i).Z),strings(currZ));
    else
        rowIndices = [];
    end
    if ~isempty(rowIndices)
        u = As(currZ);
        zEnergies = Ztable.data(i).energies(:,rowIndices); %MeV/u
        zValues = Ztable.data(i).zs(:,rowIndices);

        specEnergies = spectrum.energyBin; %MeV/u?

        dEdxFragment = dEdx(i).dEdx;
        dEdxEnergies = dEdx(i).energies; 

        dEdxInterp = interp1(dEdxEnergies,dEdxFragment,specEnergies,"linear","extrap");
        zInterp = interp1(zEnergies,zValues,specEnergies,"linear","extrap");

        alphaInterp = interp1(modelEnergies,tableEntry(i).alpha,specEnergies,"linear","extrap");
        betaInterp  = interp1(modelEnergies,tableEntry(i).beta,specEnergies,"linear","extrap");

        hC = semilogx(hDeDx,specEnergies,dEdxInterp,'LineWidth',2.5,'DisplayName',sprintf('Z = %d, A = %d',currZ,As(currZ))); hold(hDeDx,"on");
        semilogx(hZstar,specEnergies,zInterp,'LineWidth',2.5,'DisplayName',sprintf('Z = %d, A = %d',currZ,As(currZ))); hold(hZstar,"on");

        %dEdxInterp = dEdxInterp_rbetable;
        

        dividend_z_HIT = dividend_z_HIT + sum(zInterp.*dEdxInterp.*spectrum.fluenceSpectrum',2);
        denominator_z_HIT = denominator_z_HIT + sum(dEdxInterp.*spectrum.fluenceSpectrum',2);
        
        dividend_LET_HIT = dividend_LET_HIT + sum(dEdxInterp.^2.*spectrum.fluenceSpectrum',2);
        denominator_LET_HIT = denominator_LET_HIT + sum(dEdxInterp.*spectrum.fluenceSpectrum',2);
        
        dividend_alpha = dividend_alpha + sum(alphaInterp.*dEdxInterp.*spectrum.fluenceSpectrum',2);
        denominator_alpha = denominator_alpha + sum(dEdxInterp.*spectrum.fluenceSpectrum',2);

        dividend_sqrtbeta = dividend_sqrtbeta + sum(sqrt(betaInterp).*dEdxInterp.*spectrum.fluenceSpectrum',2);
        denominator_sqrtbeta = denominator_sqrtbeta + sum(dEdxInterp.*spectrum.fluenceSpectrum',2);


        %denominator = 
    end
end

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

zMix_HIT = dividend_z_HIT ./ denominator_z_HIT;
zMix_HIT(~isfinite(zMix_HIT)) = 0;
yyaxis(hMix,"left");
plot(hMix,depths_HIT,zMix_HIT,'LineWidth',2.5, 'LineStyle','--','Color','b'); hold on;
if exist('zDepthC','var')
    plot(hMix,zDepthC.depth,zDepthC.z,'LineWidth',1,'LineStyle',':','Color','k');
end
ylim(hMix,[0 25]);
ylabel(hMix,'z^* HIT');
yyaxis(hMix,"right");

zMix_MCF = dividend_z_MCF ./ denominator_z_MCF;
zMix_MCF(~isfinite(zMix_MCF)) = 0;
plot(hMix,depths_MCF,zMix_MCF,'LineWidth',1.5,'LineStyle',':');
ylabel(hMix,'z^* MCF');
xlim(hMix,[0 200]);
xlabel(hMix,'depth [mm]');
grid(hMix,"minor");
title(hMix,sprintf('Averaged z* (E_{prim} = %g MeV/u) MCF',bdEntry.energy));

hMixLET = subplot(2,2,4);
LET_HIT = dividend_LET_HIT ./ denominator_LET_HIT;
LET_HIT(~isfinite(LET_HIT)) = 0;
yyaxis(hMixLET,"left");
plot(hMixLET,depths_HIT,LET_HIT,'LineWidth',2.5,'LineStyle','--','Color','b');
ylim(hMixLET,[0 200]);
ylabel(hMixLET,'LET HIT');
yyaxis(hMixLET,"right");

LET_MCF = dividend_LET_MCF ./ denominator_LET_MCF;
LET_MCF(~isfinite(LET_MCF)) = 0;
plot(hMixLET,depths_MCF,LET_MCF,'LineWidth',1.5,'LineStyle',':');
xlim(hMixLET,[0 200]);
ylabel(hMixLET,'LET MCF');
xlabel(hMixLET,'depth [mm]');
grid(hMixLET,"minor");
title(hMixLET,sprintf('Averaged LET (E_{prim} = %g MeV/u) MCF',bdEntry.energy));

%% Plot alpha
figure;
%alpha0 = 0.0708; %0.282;
%beta0 = 0.0615;
alpha0 = tableEntry.alphaX;
beta0  = tableEntry.betaX;

alphaMKM = (alpha0 + beta0*zMix);
alphaSpec = dividend_alpha ./ denominator_alpha;

plot(depths_HIT,alphaMKM,'LineWidth',2.5); hold on;
plot(depths_HIT,alphaSpec,'LineWidth',2.5,'LineStyle',':');







