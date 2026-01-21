function [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
% matRad inverse planning wrapper function
%
% call
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln)
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   wInit:      (optional) custom weights to initialize problems
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   optimizer:  Used Optimizer Object
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

if exist('wInit','var')
    [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln,wInit);
else
    [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln);
end

%Dummy
tmpConstRBExDCheck = cellfun(@(quantity) isa(quantity,'matRad_ConstantRBExDose'), optiProb.BP.quantities, 'UniformOutput',false);
if any([tmpConstRBExDCheck{:}]) && ~isfield(dij,'RBE')
    dij.RBE = 1.1;
end

if ~isfield(pln.propOpt,'optimizer')
    %While the default optimizer is IPOPT, we can try to fallback to
    %fmincon in case it does not work for some reason
    if ~matRad_OptimizerIPOPT.IsAvailable()
        pln.propOpt.optimizer = 'fmincon';
    else
        pln.propOpt.optimizer = 'IPOPT';
    end   
end

switch pln.propOpt.optimizer
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    case 'simulannealbnd'
        optimizer = matRad_OptimizerSimulannealbnd;
    otherwise
        warning(['Optimizer ''' pln.propOpt.optimizer ''' not known! Fallback to IPOPT!']);
        optimizer = matRad_OptimizerIPOPT;
end
        
if ~optimizer.IsAvailable()
    matRad_cfg.dispError(['Optimizer ''' pln.propOpt.optimizer ''' not available!']);
end

optimizer = optimizer.optimize(wInit,optiProb,dij,cst);

wOpt = optimizer.wResult;
info = optimizer.resultInfo;

try
    resultGUI = matRad_calcCubes(wOpt,dij,[],[],cst,pln);
catch
    matRad_cfg.dispWarning('Unable to compute calcCubes');
end
resultGUI.wUnsequenced = wOpt;
resultGUI.usedOptimizer = optimizer;
resultGUI.info = info;
resultGUI.info.timePerIteration = resultGUI.info.cpu/resultGUI.info.iter;

if ~exist('computeScenarios', 'var') || isempty(computeScenarios)
    computeScenarios = 1;
end

%Robust quantities
try
    if computeScenarios
        if FLAG_ROB_OPT
            if pln.multScen.totNumScen > 1
                for i = 1:pln.multScen.totNumScen
                    scenSubIx = pln.multScen.linearMask(i,:);
                    resultGUItmp = matRad_calcCubes(wOpt,dij,pln.multScen.sub2scenIx(scenSubIx(1),scenSubIx(2),scenSubIx(3)));
                    resultGUI = matRad_appendResultGUI(resultGUI,resultGUItmp,false,sprintf('scen%d',i));
                end
            end
        end
    end
catch
    matRad_cfg.dispWarning('Unable to compute calcCubes');
end
% unblock mex files
clear mex
