%function to calculate beam optics used by mcSquare for given energy
load("carbon_HITfixedBL.mat")
energyIx = 1:255;

E = [];
d = {};
z = [];
sigmaSqIso = [];
rho = [];
sigmaT = [];

for i = energyIx
    E(i,1) = machine.data(i).energy;
    d{i,1} = machine.data(i).depths;

    focusIndex = 4;

    %calculate geometric distances and extrapolate spot size at nozzle
    SAD = machine.meta.SAD;
    z   = -(machine.data(i).initFocus.dist(focusIndex,:) - SAD);
    sigma = machine.data(i).initFocus.sigma(focusIndex,:);

    %Double Gaussian data might have a non-zero wide Gaussian,
    %adding width to the beam. We do a maximum correction here,
    %which compromises the fwhm, but seems to work better in
    %estimating the optics

    % sigmaSq_Narr = sigma.^2 + machine.data(i).sigma1(1).^2;
    % sigmaSq_Bro  = sigma.^2 + machine.data(i).sigma2(1).^2;
    % dgWeight = machine.data(i).weight(1);

    %Maximum of double gaussian
    %maxL = (1-dgWeight) ./ (2*pi*sigmaSq_Narr) + dgWeight ./ (2*pi*sigmaSq_Bro );

    %Find the sigma that corresponds to the maximum
    %sigma = sqrt(1 ./ (2*pi*maxL));

    %correct for in-air scattering with polynomial or interpolation
    % sigmaAir = polyFit(E,d);
    % sigmaAirCorrected = sqrt(sigma.^2 - sigmaAir.^2);

    % sigma_old = sigma;
    % sigma = sigmaAirCorrected;

    %square and interpolate at isocenter
    sigmaSq = sigma.^2;
    sigmaSqIso(i,1) = sqrt(interp1(z,sigmaSq,0));

    %fit Courant-Synder equation to data using ipopt, formulae
    %given in mcSquare documentation

    %fit function
    qRes = @(rho, sigmaT) (sigmaSq -  (sigmaSqIso(i,1)^2 - 2*sigmaSqIso(i,1)*rho*sigmaT.*z + sigmaT^2.*z.^2));

    % fitting for either matlab or octave
    funcs.objective = @(x) sum(qRes(x(1), x(2)).^2);
    funcs.gradient  = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso(i,1) * x(2) * z));
            2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso(i,1) * x(1) * z  - 2 * x(2) * z.^2))];

        options.lb = [-0.99, -Inf];
        options.ub = [ 0.99,  Inf];

        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs';

        %Set Default Options
        matRad_cfg = MatRad_Config.instance();
        if matRad_cfg.logLevel <= 1
            lvl = 0;
        else
            lvl = 1;
        end
        options.ipopt.print_level = lvl;

        start = [0.9; 0.1];
        [result, ~] = ipopt (start, funcs, options);
        rho(i,1)    = result(1);
        sigmaT(i,1) = result(2);

        phi{1} = @(x) sum(qRes(x(1), x(2)).^2);
        phi{2} = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso(i,1) * x(2) * z));
            2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso(i,1) * x(1) * z  - 2 * x(2) * z.^2))];

        % lb = [-0.99, -Inf];
        % ub = [ 0.99,  Inf];
        % 
        % start = [0.9; 0.1];
        % [result, ~] = sqp (start, phi, [], [], lb, ub);
        % rho    = result(1);
        % sigmaT = result(2);
end

%%
save("BeamOptics_AllFocusIx.mat","focusIndex","E","rho","sigmaSqIso","sigmaT","d")

%% Plotting

figure
subplot(3,1,1)
plot(energyIx,sigmaSqIso)
title("sigmaIso")
subplot(3,1,2)
plot(energyIx,rho)
title("rho")
subplot(3,1,3)
plot(energyIx,sigmaT)
title("sigmaT")
