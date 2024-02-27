function weightGradient = matRad_objectiveGradient(optiProb,w,dij,cst)
% matRad IPOPT callback: gradient function for inverse planning
% supporting mean dose objectives, EUD objectives, squared overdosage,
% squared underdosage, squared deviation and DVH objectives
%
% call
%   g = matRad_gradFuncWrapper(optiProb,w,dij,cst)
%
% input
%   optiProb: option struct defining the type of optimization
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%
% output
%   g: gradient of objective function
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();
% DOSE PROJECTION
bxidx = 1; %modality  bixel index

% Obtain cumulative dose cube 
[d{1}]    = deal(zeros(dij.doseGrid.numOfVoxels,1));
[dD{1}]   = deal(zeros(dij.doseGrid.numOfVoxels,1));
[dLD{1}]   = deal(zeros(dij.doseGrid.numOfVoxels,1));

for mod = 1: length(dij.original_Dijs)
    wt = [];
    % split the w for current modality
    STrepmat = (~dij.spatioTemp(mod) + dij.spatioTemp(mod)*dij.numOfSTscen(mod));
    wt = reshape(w(bxidx: bxidx+STrepmat*dij.original_Dijs{mod}.totalNumOfBixels-1),[dij.original_Dijs{mod}.totalNumOfBixels,STrepmat]);
    % get current dose / effect / RBExDose vector
    optiProb.BP.compute(dij.original_Dijs{mod},wt);
    dTmp = optiProb.BP.GetResult();

    % also get probabilistic quantities (nearly no overhead if empty)
    [dExp,dOmega] = optiProb.BP.GetResultProb();                        % NOTE: not sure what exactly to do for the dOmegas 

    % Manage other quantites
    if ~isempty(optiProb.dirtyDoseBP)
        optiProb.dirtyDoseBP.compute(dij.original_Dijs{mod},wt);
        dDTmp = optiProb.dirtyDoseBP.GetResult();

        % also get probabilistic quantities (nearly no overhead if empty)
        [dDExp,dDOmega] = optiProb.dirtyDoseBP.GetResultProb();
    end

    if ~isempty(optiProb.mLETDoseBP)
        optiProb.mLETDoseBP.compute(dij.original_Dijs{mod},wt);
        dLDTmp = optiProb.mLETDoseBP.GetResult();

        % also get probabilistic quantities (nearly no overhead if empty)
        [dLDExp,dLDOmega] = optiProb.mLETDoseBP.GetResultProb();

    end

    %Index of nxt modality
    bxidx = bxidx + STrepmat*dij.original_Dijs{mod}.totalNumOfBixels;
    % Accumulat Dose for all scenarios FOR FUTURE REVIEW ON HOW TO COMBINE
    % DIFFFERENT UNCERTAINTY SCENARIOS FOR DIFFERENT MODALITIES
    % currently for ST optimization 
    for scen = 1:numel(dTmp)
        d{scen} = d{scen} + sum(dTmp{scen}.*dij.STfractions{mod}',2);
    end
    if ~isempty(optiProb.dirtyDoseBP)
        for scen = 1:numel(dDTmp)
            dD{scen} = dD{scen} + sum(dDTmp{scen}.*dij.STfractions{mod}',2);
        end
    end
    if ~isempty(optiProb.mLETDoseBP)
        for scen = 1:numel(dLDTmp)
            dLD{scen} = dLD{scen} + sum(dLDTmp{scen}.*dij.STfractions{mod}',2);
        end
    end
end

    % get the used scenarios
    useScen  = optiProb.BP.scenarios;
    scenProb = optiProb.BP.scenarioProb;

    if ~isempty(optiProb.dirtyDoseBP)
        dirtyDoseGradient          = cell(size(dij.original_Dijs{mod}.dirtyDose));
        dirtyDoseGradient(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
    end
    if ~isempty(optiProb.mLETDoseBP)
        mLETDoseGradient          = cell(size(dij.original_Dijs{mod}.mLETDose));
        mLETDoseGradient(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
    end

    % retrieve matching 4D scenarios
    fullScen      = cell(ndims(dTmp),1);
    [fullScen{:}] = ind2sub(size(dTmp),useScen);
    contourScen   = fullScen{1};

    doseGradient          = cell(size(dij.original_Dijs{mod}.physicalDose));
    doseGradient(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};

    %For probabilistic optimization
    vOmega = 0;

    %For COWC
    f_COWC = zeros(size(dij.original_Dijs{mod}.physicalDose));

    % compute objective function for every VOI.
    for  i = 1:size(cst,1)

        % Only take OAR or target VOI.
        if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

            % loop over the number of constraints and objectives for the current VOI
            for j = 1:numel(cst{i,6})

                %Get current optimization function
                objective = cst{i,6}{j};

                % only perform gradient computations for objectives
                if isa(objective,'DoseObjectives.matRad_DoseObjective')

                    % rescale dose parameters to biological optimization quantity if required

                    % Obtain the fraction-size dose prescription to be
                    % converted to biological prescription. When performing
                    % effect-based optimization, this step ensures that the
                    % dose-effect conversion is always performed on a
                    % fraction-size dose prescription and only after scaled back to
                    % total plan effect.

                    % First scale the dose prescription to fraction size
                    doseParameter = objective.getDoseParameters();
                    objective = objective.setDoseParameters(doseParameter./dij.totalNumOfFractions);

                    % Compute the biological parameters. This function is only effective when using effect-based
                    % optimization
 
                    objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
               
                    % Scale back the biological dose prescription to total plan
                    % size
                    doseParameter = objective.getDoseParameters();
                    objective = objective.setDoseParameters(doseParameter.*dij.totalNumOfFractions); 

                    % retrieve the robustness type
                    robustness = objective.robustness;


                    switch robustness
                        case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                d_i = d{ixScen}(cst{i,4}{ixContour});
                                %add to dose gradient
                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + objective.penalty * objective.computeDoseObjectiveGradient(d_i);
                            end
                        case 'STOCH' % perform stochastic optimization with weighted / random scenarios
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = d{ixScen}(cst{i,4}{ixContour});

                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + ...
                                    (objective.penalty * objective.computeDoseObjectiveGradient(d_i) * scenProb(s));

                            end

                        case 'PROB' % use the expectation value and the integral variance influence matrix
                            %First check the speficic cache for probabilistic
                            if ~exist('doseGradientExp','var')
                                doseGradientExp{1} = zeros(dij.doseGrid.numOfVoxels,1);
                            end

                            d_i = dExp{1}(cst{i,4}{1});

                            doseGradientExp{1}(cst{i,4}{1}) = doseGradientExp{1}(cst{i,4}{1}) + objective.penalty * objective.computeDoseObjectiveGradient(d_i);

                            p = objective.penalty/numel(cst{i,4}{1});

                            vOmega = vOmega + p * dOmega{i,1};

                        case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                            contourIx = unique(contourScen);
                            if ~isscalar(contourIx)
                                % voxels need to be tracked through the 4D CT,
                                % not yet implemented
                                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                            end

                            % prepare min/max dose vector for voxel-wise worst case
                            if ~exist('d_tmp','var')
                                d_tmp = [d{useScen}];
                            end

                            d_Scen = d_tmp(cst{i,4}{contourIx},:);
                            [d_max,max_ix] = max(d_Scen,[],2);
                            [d_min,min_ix] = min(d_Scen,[],2);

                            if isequal(cst{i,3},'OAR')
                                d_i = d_max;
                            elseif isequal(cst{i,3},'TARGET')
                                d_i = d_min;
                            end

                            if any(isnan(d_i))
                                matRad_cfg.dispWarning('%d NaN values in gradient.',numel(isnan(d_i)));
                            end

                            deltaTmp = objective.penalty * objective.computeDoseObjectiveGradient(d_i);

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                if isequal(cst{i,3},'OAR')
                                    currWcIx = double(max_ix == s);
                                elseif isequal(cst{i,3},'TARGET')
                                    currWcIx = double(min_ix == s);
                                end

                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                            end

                        case 'VWWC_INV'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                            contourIx = unique(contourScen);
                            if ~isscalar(contourIx)
                                % voxels need to be tracked through the 4D CT,
                                % not yet implemented
                                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                            end

                            % prepare min/max dose vector for voxel-wise worst case
                            if ~exist('d_tmp','var')
                                d_tmp = [d{useScen}];
                            end

                            d_Scen = d_tmp(cst{i,4}{1},:);
                            [d_max,max_ix] = max(d_Scen,[],2);
                            [d_min,min_ix] = min(d_Scen,[],2);

                            if isequal(cst{i,3},'OAR')
                                d_i = d_min;
                            elseif isequal(cst{i,3},'TARGET')
                                d_i = d_max;
                            end

                            if any(isnan(d_i))
                                matRad_cfg.dispWarning('%d NaN values in gradFuncWrapper.',numel(isnan(d_i)));
                            end

                            deltaTmp = objective.penalty * objective.computeDoseObjectiveGradient(d_i);

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                if isequal(cst{i,3},'OAR')
                                    currWcIx = double(min_ix == s);
                                elseif isequal(cst{i,3},'TARGET')
                                    currWcIx = double(max_ix == s);
                                end

                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                            end

                        case 'COWC' % composite worst case consideres ovarall the worst objective function value
                            %First check the speficic cache for COWC
                            if ~exist('delta_COWC','var')
                                delta_COWC         = cell(size(doseGradient));
                                delta_COWC(useScen)    = {zeros(dij.doseGrid.numOfVoxels,1)};
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = d{ixScen}(cst{i,4}{ixContour});

                                f_COWC(ixScen) = f_COWC(ixScen) + objective.penalty * objective.computeDoseObjectiveFunction(d_i);
                                delta_COWC{ixScen}(cst{i,4}{ixContour}) = delta_COWC{ixScen}(cst{i,4}{ixContour}) + objective.penalty * objective.computeDoseObjectiveGradient(d_i);
                            end

                        case 'OWC' % objective-wise worst case consideres the worst individual objective function value
                            %First check the speficic cache for COWC
                            f_OWC = zeros(size(doseGradient));

                            if ~exist('delta_OWC','var')
                                delta_OWC = cell(size(doseGradient));
                                delta_OWC(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = d{ixScen}(cst{i,4}{ixContour});

                                f_OWC(ixScen) = objective.penalty * objective.computeDoseObjectiveFunction(d_i);

                                delta_OWC{ixScen}(cst{i,4}{ixContour}) = objective.penalty * objective.computeDoseObjectiveGradient(d_i);

                            end

                            switch optiProb.useMaxApprox
                                case 'logsumexp'
                                    [~,fGrad] = optiProb.logSumExp(f_OWC);
                                case 'pnorm'
                                    [~,fGrad] = optiProb.pNorm(f_OWC,numel(useScen));
                                case 'none'
                                    [~,ix] = max(f_OWC(:));
                                    fGrad = zeros(size(f_OWC));
                                    fGrad(ix) = 1;
                                case 'otherwise'
                                    matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                                    [~,ix] = max(f_OWC(:));
                                    fGrad = zeros(size(f_OWC));
                                    fGrad(ix) = 1;
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                if fGrad(ixScen ) ~= 0
                                    doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + fGrad(ixScen)*delta_OWC{ixScen}(cst{i,4}{ixContour});
                                end
                            end

                        otherwise
                            matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);

                    end  %robustness type
                end  % objective check

                % if ~isempty(optiProb.dirtyDoseBP)
                %     rowCst = size(cst);
                %     for a = 1 : rowCst(1,1)
                %         num = size(cst{a,6});
                %         for b = 1 : num(1,2)
                %             if isa(cst{a,6}{b},'DirtyDoseObjectives.matRad_DirtyDoseObjective')
                %                 isa(objective,'DirtyDoseObjectives.matRad_DirtyDoseObjective')
                %             end
                %         end
                %     end
                % end

                if ~isempty(optiProb.dirtyDoseBP) && isa(objective,'DirtyDoseObjectives.matRad_DirtyDoseObjective')

                    % rescale dose parameters to biological optimization quantity if required

                    % Obtain the fraction-size dose prescription to be
                    % converted to biological prescription. When performing
                    % effect-based optimization, this step ensures that the
                    % dose-effect conversion is always performed on a
                    % fraction-size dose prescription and only after scaled back to
                    % total plan effect.

                    % First scale the dose prescription to fraction size
                    doseParameter = objective.getDoseParameters();
                    objective = objective.setDoseParameters(doseParameter./dij.totalNumOfFractions);

                    % Compute the biological parameters. This function is only effective when using effect-based
                    % optimization
 
                    %objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
               
                    % Scale back the biological dose prescription to total plan
                    % size
                    doseParameter = objective.getDoseParameters();
                    objective = objective.setDoseParameters(doseParameter.*dij.totalNumOfFractions); 

                    % retrieve the robustness type
                    robustness = objective.robustness;


                    switch robustness
                        case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                d_i = dD{ixScen}(cst{i,4}{ixContour});
                                %add to dose gradient
                                dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) = dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) + objective.penalty * objective.computeDirtyDoseObjectiveGradient(d_i);
                            end
                        case 'STOCH' % perform stochastic optimization with weighted / random scenarios
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = dD{ixScen}(cst{i,4}{ixContour});

                                dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) = dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) + ...
                                    (objective.penalty * objective.computeDirtyDoseObjectiveGradient(d_i) * scenProb(s));

                            end

                        case 'PROB' % use the expectation value and the integral variance influence matrix
                            %First check the speficic cache for probabilistic
                            if ~exist('dirtyDoseGradientExp','var')
                                dirtyDoseGradientExp{1} = zeros(dij.doseGrid.numOfVoxels,1);
                            end

                            d_i = dDExp{1}(cst{i,4}{1});

                            dirtyDoseGradientExp{1}(cst{i,4}{1}) = dirtyDoseGradientExp{1}(cst{i,4}{1}) + objective.penalty * objective.computeDirtyDoseObjectiveGradient(d_i);

                            p = objective.penalty/numel(cst{i,4}{1});

                            vOmega = vOmega + p * dDOmega{i,1};

                        case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                            contourIx = unique(contourScen);
                            if ~isscalar(contourIx)
                                % voxels need to be tracked through the 4D CT,
                                % not yet implemented
                                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                            end

                            % prepare min/max dose vector for voxel-wise worst case
                            if ~exist('d_tmp','var')
                                d_tmp = [dD{useScen}];
                            end

                            d_Scen = d_tmp(cst{i,4}{contourIx},:);
                            [d_max,max_ix] = max(d_Scen,[],2);
                            [d_min,min_ix] = min(d_Scen,[],2);

                            if isequal(cst{i,3},'OAR')
                                d_i = d_max;
                            elseif isequal(cst{i,3},'TARGET')
                                d_i = d_min;
                            end

                            if any(isnan(d_i))
                                matRad_cfg.dispWarning('%d NaN values in gradient.',numel(isnan(d_i)));
                            end

                            deltaTmp = objective.penalty * objective.computeDirtyDoseObjectiveGradient(d_i);

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                if isequal(cst{i,3},'OAR')
                                    currWcIx = double(max_ix == s);
                                elseif isequal(cst{i,3},'TARGET')
                                    currWcIx = double(min_ix == s);
                                end

                                dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) = dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                            end

                        case 'VWWC_INV'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                            contourIx = unique(contourScen);
                            if ~isscalar(contourIx)
                                % voxels need to be tracked through the 4D CT,
                                % not yet implemented
                                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                            end

                            % prepare min/max dose vector for voxel-wise worst case
                            if ~exist('d_tmp','var')
                                d_tmp = [dD{useScen}];
                            end

                            d_Scen = d_tmp(cst{i,4}{1},:);
                            [d_max,max_ix] = max(d_Scen,[],2);
                            [d_min,min_ix] = min(d_Scen,[],2);

                            if isequal(cst{i,3},'OAR')
                                d_i = d_min;
                            elseif isequal(cst{i,3},'TARGET')
                                d_i = d_max;
                            end

                            if any(isnan(d_i))
                                matRad_cfg.dispWarning('%d NaN values in gradFuncWrapper.',numel(isnan(d_i)));
                            end

                            deltaTmp = objective.penalty * objective.computeDirtyDoseObjectiveGradient(d_i);

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                if isequal(cst{i,3},'OAR')
                                    currWcIx = double(min_ix == s);
                                elseif isequal(cst{i,3},'TARGET')
                                    currWcIx = double(max_ix == s);
                                end

                                dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) = dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                            end

                        case 'COWC' % composite worst case consideres ovarall the worst objective function value
                            %First check the speficic cache for COWC
                            if ~exist('delta_COWC','var')
                                delta_COWC         = cell(size(dirtyDoseGradient));
                                delta_COWC(useScen)    = {zeros(dij.doseGrid.numOfVoxels,1)};
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = dD{ixScen}(cst{i,4}{ixContour});

                                f_COWC(ixScen) = f_COWC(ixScen) + objective.penalty * objective.computeDirtyDoseObjectiveFunction(d_i);
                                delta_COWC{ixScen}(cst{i,4}{ixContour}) = delta_COWC{ixScen}(cst{i,4}{ixContour}) + objective.penalty * objective.computeDirtyDoseObjectiveGradient(d_i);
                            end

                        case 'OWC' % objective-wise worst case consideres the worst individual objective function value
                            %First check the speficic cache for COWC
                            f_OWC = zeros(size(dirtyDoseGradient));

                            if ~exist('delta_OWC','var')
                                delta_OWC = cell(size(dirtyDoseGradient));
                                delta_OWC(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = dD{ixScen}(cst{i,4}{ixContour});

                                f_OWC(ixScen) = objective.penalty * objective.computeDirtyDoseObjectiveFunction(d_i);

                                delta_OWC{ixScen}(cst{i,4}{ixContour}) = objective.penalty * objective.computeDirtyDoseObjectiveGradient(d_i);

                            end

                            switch optiProb.useMaxApprox
                                case 'logsumexp'
                                    [~,fGrad] = optiProb.logSumExp(f_OWC);
                                case 'pnorm'
                                    [~,fGrad] = optiProb.pNorm(f_OWC,numel(useScen));
                                case 'none'
                                    [~,ix] = max(f_OWC(:));
                                    fGrad = zeros(size(f_OWC));
                                    fGrad(ix) = 1;
                                case 'otherwise'
                                    matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                                    [~,ix] = max(f_OWC(:));
                                    fGrad = zeros(size(f_OWC));
                                    fGrad(ix) = 1;
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                if fGrad(ixScen ) ~= 0
                                    dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) = dirtyDoseGradient{ixScen}(cst{i,4}{ixContour}) + fGrad(ixScen)*delta_OWC{ixScen}(cst{i,4}{ixContour});
                                end
                            end

                        otherwise
                            matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);

                    end  %robustness type
                end  % objective check

                % if ~isempty(optiProb.mLETDoseBP)
                %     rowCst = size(cst);
                %     for a = 1 : rowCst(1,1)
                %         num = size(cst{a,6});
                %         for b = 1 : num(1,2)
                %             if isa(cst{a,6}{b},'mLETDoseObjectives.matRad_mLETDoseObjective')
                %                 isa(objective,'mLETDoseObjectives.matRad_mLETDoseObjective')
                %             end
                %         end
                %     end
                % end

                if ~isempty(optiProb.mLETDoseBP) && isa(objective,'mLETDoseObjectives.matRad_mLETDoseObjective')

                    % rescale dose parameters to biological optimization quantity if required

                    % Obtain the fraction-size dose prescription to be
                    % converted to biological prescription. When performing
                    % effect-based optimization, this step ensures that the
                    % dose-effect conversion is always performed on a
                    % fraction-size dose prescription and only after scaled back to
                    % total plan effect.

                    % First scale the dose prescription to fraction size
                    doseParameter = objective.getDoseParameters();
                    objective = objective.setDoseParameters(doseParameter./dij.totalNumOfFractions);

                    % Compute the biological parameters. This function is only effective when using effect-based
                    % optimization
 
                    %objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
               
                    % Scale back the biological dose prescription to total plan
                    % size
                    doseParameter = objective.getDoseParameters();
                    objective = objective.setDoseParameters(doseParameter.*dij.totalNumOfFractions); 

                    % retrieve the robustness type
                    robustness = objective.robustness;


                    switch robustness
                        case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                d_i = dLD{ixScen}(cst{i,4}{ixContour});
                                %add to dose gradient
                                mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) = mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) + objective.penalty * objective.computemLETDoseObjectiveGradient(d_i);
                            end
                        case 'STOCH' % perform stochastic optimization with weighted / random scenarios
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = dLD{ixScen}(cst{i,4}{ixContour});

                                mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) = mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) + ...
                                    (objective.penalty * objective.computemLETDoseObjectiveGradient(d_i) * scenProb(s));

                            end

                        case 'PROB' % use the expectation value and the integral variance influence matrix
                            %First check the speficic cache for probabilistic
                            if ~exist('mLETDoseGradientExp','var')
                                mLETDoseGradientExp{1} = zeros(dij.doseGrid.numOfVoxels,1);
                            end

                            d_i = dLDExp{1}(cst{i,4}{1});

                            mLETDoseGradientExp{1}(cst{i,4}{1}) = mLETDoseGradientExp{1}(cst{i,4}{1}) + objective.penalty * objective.computemLETDoseObjectiveGradient(d_i);

                            p = objective.penalty/numel(cst{i,4}{1});

                            vOmega = vOmega + p * dLDOmega{i,1};

                        case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                            contourIx = unique(contourScen);
                            if ~isscalar(contourIx)
                                % voxels need to be tracked through the 4D CT,
                                % not yet implemented
                                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                            end

                            % prepare min/max dose vector for voxel-wise worst case
                            if ~exist('d_tmp','var')
                                d_tmp = [dLD{useScen}];
                            end

                            d_Scen = d_tmp(cst{i,4}{contourIx},:);
                            [d_max,max_ix] = max(d_Scen,[],2);
                            [d_min,min_ix] = min(d_Scen,[],2);

                            if isequal(cst{i,3},'OAR')
                                d_i = d_max;
                            elseif isequal(cst{i,3},'TARGET')
                                d_i = d_min;
                            end

                            if any(isnan(d_i))
                                matRad_cfg.dispWarning('%d NaN values in gradient.',numel(isnan(d_i)));
                            end

                            deltaTmp = objective.penalty * objective.computemLETDoseObjectiveGradient(d_i);

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                if isequal(cst{i,3},'OAR')
                                    currWcIx = double(max_ix == s);
                                elseif isequal(cst{i,3},'TARGET')
                                    currWcIx = double(min_ix == s);
                                end

                                mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) = mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                            end

                        case 'VWWC_INV'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                            contourIx = unique(contourScen);
                            if ~isscalar(contourIx)
                                % voxels need to be tracked through the 4D CT,
                                % not yet implemented
                                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                            end

                            % prepare min/max dose vector for voxel-wise worst case
                            if ~exist('d_tmp','var')
                                d_tmp = [dLD{useScen}];
                            end

                            d_Scen = d_tmp(cst{i,4}{1},:);
                            [d_max,max_ix] = max(d_Scen,[],2);
                            [d_min,min_ix] = min(d_Scen,[],2);

                            if isequal(cst{i,3},'OAR')
                                d_i = d_min;
                            elseif isequal(cst{i,3},'TARGET')
                                d_i = d_max;
                            end

                            if any(isnan(d_i))
                                matRad_cfg.dispWarning('%d NaN values in gradFuncWrapper.',numel(isnan(d_i)));
                            end

                            deltaTmp = objective.penalty * objective.computemLETDoseObjectiveGradient(d_i);

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                if isequal(cst{i,3},'OAR')
                                    currWcIx = double(min_ix == s);
                                elseif isequal(cst{i,3},'TARGET')
                                    currWcIx = double(max_ix == s);
                                end

                                mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) = mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                            end

                        case 'COWC' % composite worst case consideres ovarall the worst objective function value
                            %First check the speficic cache for COWC
                            if ~exist('delta_COWC','var')
                                delta_COWC         = cell(size(mLETDoseGradient));
                                delta_COWC(useScen)    = {zeros(dij.doseGrid.numOfVoxels,1)};
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = dLD{ixScen}(cst{i,4}{ixContour});

                                f_COWC(ixScen) = f_COWC(ixScen) + objective.penalty * objective.computemLETDoseObjectiveFunction(d_i);
                                delta_COWC{ixScen}(cst{i,4}{ixContour}) = delta_COWC{ixScen}(cst{i,4}{ixContour}) + objective.penalty * objective.computemLETDoseObjectiveGradient(d_i);
                            end

                        case 'OWC' % objective-wise worst case consideres the worst individual objective function value
                            %First check the speficic cache for COWC
                            f_OWC = zeros(size(mLETDoseGradient));

                            if ~exist('delta_OWC','var')
                                delta_OWC = cell(size(mLETDoseGradient));
                                delta_OWC(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);

                                d_i = d{ixScen}(cst{i,4}{ixContour});

                                f_OWC(ixScen) = objective.penalty * objective.computemLETDoseObjectiveFunction(d_i);

                                delta_OWC{ixScen}(cst{i,4}{ixContour}) = objective.penalty * objective.computemLETDoseObjectiveGradient(d_i);

                            end

                            switch optiProb.useMaxApprox
                                case 'logsumexp'
                                    [~,fGrad] = optiProb.logSumExp(f_OWC);
                                case 'pnorm'
                                    [~,fGrad] = optiProb.pNorm(f_OWC,numel(useScen));
                                case 'none'
                                    [~,ix] = max(f_OWC(:));
                                    fGrad = zeros(size(f_OWC));
                                    fGrad(ix) = 1;
                                case 'otherwise'
                                    matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                                    [~,ix] = max(f_OWC(:));
                                    fGrad = zeros(size(f_OWC));
                                    fGrad(ix) = 1;
                            end

                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                if fGrad(ixScen ) ~= 0
                                    mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) = mLETDoseGradient{ixScen}(cst{i,4}{ixContour}) + fGrad(ixScen)*delta_OWC{ixScen}(cst{i,4}{ixContour});
                                end
                            end

                        otherwise
                            matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);

                    end  %robustness type
                end  % objective check

            end %objective loop
        end %empty check
    end %cst structure loop

    if exist('delta_COWC','var')
        switch optiProb.useMaxApprox
            case 'logsumexp'
                [~,fGrad] = optiProb.logSumExp(f_COWC);
            case 'pnorm'
                [~,fGrad] = optiProb.pNorm(f_COWC,numel(useScen));
            case 'none'
                [~,ixCurrWC] = max(f_COWC(:));
                fGrad = zeros(size(f_COWC));
                fGrad(ixCurrWC) = 1;
            case 'otherwise'
                matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                [~,ixCurrWC] = max(f_COWC(:));
                fGrad = zeros(size(f_COWC));
                fGrad(ixCurrWC) = 1;
        end

        for s = 1:numel(useScen)
            ixScen = useScen(s);
            if fGrad(ixScen) ~= 0
                if ~isempty(optiProb.BP)
                    doseGradient{ixScen} = doseGradient{ixScen} + fGrad(ixScen)*delta_COWC{ixScen};
                end
                if ~isempty(optiProb.dirtyDoseBP)
                    dirtyDoseGradient{ixScen} = dirtyDoseGradient{ixScen} + fGrad{ixScen} * delta_COWC{ixScen};
                end
                if ~isempty(optiProb.mLETDoseBP)
                    mLETDoseGradient{ixScen} = mLETDoseGradient{ixScen} + fGrad{ixScen} * delta_COWC{ixScen};
                end
            end
        end

    end

bxidx = 1;
g = cell(numel(useScen),1);
if ~isempty(optiProb.dirtyDoseBP)
    gd = cell(numel(useScen),1);
end
if ~isempty(optiProb.mLETDoseBP)
    gL = cell(numel(useScen),1);
end

for mod = 1: length(dij.original_Dijs)
    wt = [];

    % split the w and g for current modality
    STrepmat = (~dij.spatioTemp(mod) + dij.spatioTemp(mod)*dij.numOfSTscen(mod));
    wt = reshape(w(bxidx: bxidx+STrepmat*dij.original_Dijs{mod}.totalNumOfBixels-1),[dij.original_Dijs{mod}.totalNumOfBixels,STrepmat]);
    
    optiProb.BP.computeGradient(dij.original_Dijs{mod},doseGradient,wt);
    gt = optiProb.BP.GetGradient();
                     % review for ST optimization 
     for s = 1:numel(useScen)
         gt{s} = gt{s}*dij.STfractions{mod};
         g{s} = [g{s}; gt{s}];
     end
     bxidx = bxidx + STrepmat*dij.original_Dijs{mod}.totalNumOfBixels;

     if ~isempty(optiProb.dirtyDoseBP)
         optiProb.dirtyDoseBP.computeGradient(dij.original_Dijs{mod},dirtyDoseGradient,wt);
         gDD = optiProb.dirtyDoseBP.GetGradient();

         for s = 1:numel(useScen)
             gDD{s} = gDD{s}*dij.STfractions{mod};
             gd{s} = [gd{s}; gDD{s}];
         end
         
     end

    if ~isempty(optiProb.mLETDoseBP) 
        optiProb.mLETDoseBP.computeGradient(dij.original_Dijs{mod},mLETDoseGradient,wt);
        gLD = optiProb.mLETDoseBP.GetGradient();
        
        for s = 1:numel(useScen)
            gLD{s} = gLD{s}*dij.STfractions{mod};
            gL{s} = [gL{s}; gLD{s}];
        end
           
    end

end

weightGradient = zeros(dij.totalNumOfBixels,1);
for s = 1:numel(useScen)
    weightGradient = weightGradient + g{useScen(s)};
end

if vOmega ~= 0
    optiProb.BP.computeGradientProb(dij.original_Dijs{mod},doseGradientExp,vOmega,w);
    gProb = optiProb.BP.GetGradientProb();

    %Only implemented for first scenario now
    weightGradient = weightGradient + gProb{1};
end
if ~isempty(optiProb.dirtyDoseBP)
   
    % weightGradient = zeros(dij.totalNumOfBixels,1);
    for s = 1:numel(useScen)
        weightGradient = weightGradient + gd{useScen(s)};
    end

    if vOmega ~= 0
        optiProb.dirtyDoseBP.computeGradientProb(dij.original_Dijs{mod},dirtyDoseGradientExp,vOmega,w);
        gdProb = optiProb.dirtyDoseBP.GetGradientProb();

        %Only implemented for first scenario now
        weightGradient = weightGradient + gdProb{1};
    end
end
if ~isempty(optiProb.mLETDoseBP)
    % weightGradient = zeros(dij.totalNumOfBixels,1);
    for s = 1:numel(useScen)
        weightGradient = weightGradient + gL{useScen(s)};
    end

    if vOmega ~= 0
        optiProb.mLETDoseBP.computeGradientProb(dij.original_Dijs{mod},mLETDoseGradientExp,vOmega,w);
        gLProb = optiProb.mLETDoseBP.GetGradientProb();

        %Only implemented for first scenario now
        weightGradient = weightGradient + gLProb{1};
    end
end

% code snippet to check the gradient
    gradientChecker = 0;
if gradientChecker == 1
    f =  matRad_objectiveFunction(optiProb,w,dij,cst);
    epsilon = 1e-6;
    ix = unique(randi([1 numel(w)],1,5));
    
    for i=ix
        wInit = w;
        wInit(i) = wInit(i) + epsilon;
        fDel= matRad_objectiveFunction(optiProb,wInit,dij,cst);
        numGrad = (fDel - f)/epsilon;
        diff = (numGrad/weightGradient(i) - 1)*100;
        fprintf(['grad val #' num2str(i) ' - rel diff numerical and analytical gradient = ' num2str(diff) '\n']);
        %fprintf([' any nan or zero for photons' num2str(sum(isnan(glog{1}))) ',' num2str(sum(~logical(glog{1}))) ' for protons: ' num2str(sum(isnan(glog{2}))) ',' num2str(sum(~logical(glog{2}))) '\n']);
    end
    
end
end
