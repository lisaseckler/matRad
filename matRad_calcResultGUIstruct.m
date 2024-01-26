function resultGUI = matRad_calcResultGUIstruct(resultGUI)
% creates resultGUI struct by adding the fields in the resultGUI cell

resultGUI_cell = resultGUI;
resultGUI = struct();

if isfield(resultGUI_cell{1,1},'LET') && isfield(resultGUI_cell{1,2},'LET')
    resultGUI.LET = resultGUI_cell{1,1}.LET + resultGUI_cell{1,2}.LET;
elseif isfield(resultGUI_cell{1,1},'LET') && ~isfield(resultGUI_cell{1,2},'LET')
    resultGUI.LET1 = resultGUI_cell{1,1}.LET;
elseif ~isfield(resultGUI_cell{1,1},'LET') && isfield(resultGUI_cell{1,2},'LET')
    resultGUI.LET2 = resultGUI_cell{1,2}.LET;
end

if isfield(resultGUI_cell{1,1},'RBE') && isfield(resultGUI_cell{1,2},'RBE')
    resultGUI.RBE = resultGUI_cell{1,1}.RBE + resultGUI_cell{1,2}.RBE;
elseif isfield(resultGUI_cell{1,1},'RBE') && ~isfield(resultGUI_cell{1,2},'RBE')
    resultGUI.RBE1 = resultGUI_cell{1,1}.RBE;
elseif ~isfield(resultGUI_cell{1,1},'RBE') && isfield(resultGUI_cell{1,2},'RBE')
    resultGUI.RBE2 = resultGUI_cell{1,2}.RBE;
end

if isfield(resultGUI_cell{1,1},'RBExD') && isfield(resultGUI_cell{1,2},'RBExD')
    resultGUI.RBExD = resultGUI_cell{1,1}.RBExD + resultGUI_cell{1,2}.RBExD;
elseif isfield(resultGUI_cell{1,1},'RBExD') && ~isfield(resultGUI_cell{1,2},'RBExD')
    resultGUI.RBExD1 = resultGUI_cell{1,1}.RBExD;
elseif ~isfield(resultGUI_cell{1,1},'RBExD') && isfield(resultGUI_cell{1,2},'RBExD')
    resultGUI.RBExD2 = resultGUI_cell{1,2}.RBExD;
end

if isfield(resultGUI_cell{1,1},'SqrtBetaDoseCube') && isfield(resultGUI_cell{1,2},'SqrtBetaDoseCube')
    resultGUI.SqrtBetaDoseCube = resultGUI_cell{1,1}.SqrtBetaDoseCube + resultGUI_cell{1,2}.SqrtBetaDoseCube;
elseif isfield(resultGUI_cell{1,1},'SqrtBetaDoseCube') && ~isfield(resultGUI_cell{1,2},'SqrtBetaDoseCube')
    resultGUI.SqrtBetaDoseCube1 = resultGUI_cell{1,1}.SqrtBetaDoseCube;
elseif ~isfield(resultGUI_cell{1,1},'SqrtBetaDoseCube') && isfield(resultGUI_cell{1,2},'SqrtBetaDoseCube')
    resultGUI.SqrtBetaDoseCube2 = resultGUI_cell{1,2}.SqrtBetaDoseCube;
end

if isfield(resultGUI_cell{1,1},'alpha') && isfield(resultGUI_cell{1,2},'alpha')
    resultGUI.alpha = resultGUI_cell{1,1}.alpha + resultGUI_cell{1,2}.alpha;
elseif isfield(resultGUI_cell{1,1},'alpha') && ~isfield(resultGUI_cell{1,2},'alpha')
    resultGUI.alpha1 = resultGUI_cell{1,1}.alpha;
elseif ~isfield(resultGUI_cell{1,1},'alpha') && isfield(resultGUI_cell{1,2},'alpha')
    resultGUI.alpha2 = resultGUI_cell{1,2}.alpha;
end

if isfield(resultGUI_cell{1,1},'alphaDoseCube') && isfield(resultGUI_cell{1,2},'alphaDoseCube')
    resultGUI.alphaDoseCube = resultGUI_cell{1,1}.alphaDoseCube + resultGUI_cell{1,2}.alphaDoseCube;
elseif isfield(resultGUI_cell{1,1},'alphaDoseCube') && ~isfield(resultGUI_cell{1,2},'alphaDoseCube')
    resultGUI.alphaDoseCube1 = resultGUI_cell{1,1}.alphaDoseCube;
elseif ~isfield(resultGUI_cell{1,1},'alphaDoseCube') && isfield(resultGUI_cell{1,2},'alphaDoseCube')
    resultGUI.alphaDoseCube2 = resultGUI_cell{1,2}.alphaDoseCube;
end

if isfield(resultGUI_cell{1,1},'beta') && isfield(resultGUI_cell{1,2},'beta')
    resultGUI.beta = resultGUI_cell{1,1}.beta + resultGUI_cell{1,2}.beta;
elseif isfield(resultGUI_cell{1,1},'beta') && ~isfield(resultGUI_cell{1,2},'beta')
    resultGUI.beta1 = resultGUI_cell{1,1}.beta;
elseif ~isfield(resultGUI_cell{1,1},'beta') && isfield(resultGUI_cell{1,2},'beta')
    resultGUI.beta2 = resultGUI_cell{1,2}.beta;
end

if isfield(resultGUI_cell{1,1},'cleanDose') && isfield(resultGUI_cell{1,2},'cleanDose')
    resultGUI.cleanDose = resultGUI_cell{1,1}.cleanDose + resultGUI_cell{1,2}.cleanDose;
elseif isfield(resultGUI_cell{1,1},'cleanDose') && ~isfield(resultGUI_cell{1,2},'cleanDose')
    resultGUI.cleanDose1 = resultGUI_cell{1,1}.cleanDose;
elseif ~isfield(resultGUI_cell{1,1},'cleanDose') && isfield(resultGUI_cell{1,2},'cleanDose')
    resultGUI.cleanDose2 = resultGUI_cell{1,2}.cleanDose;
end

if isfield(resultGUI_cell{1,1},'dirtyDose') && isfield(resultGUI_cell{1,2},'dirtyDose')
    resultGUI.dirtyDose = resultGUI_cell{1,1}.dirtyDose + resultGUI_cell{1,2}.dirtyDose;
elseif isfield(resultGUI_cell{1,1},'dirtyDose') && ~isfield(resultGUI_cell{1,2},'dirtyDose')
    resultGUI.dirtyDose1 = resultGUI_cell{1,1}.dirtyDose;
elseif ~isfield(resultGUI_cell{1,1},'dirtyDose') && isfield(resultGUI_cell{1,2},'dirtyDose')
    resultGUI.dirtyDose2 = resultGUI_cell{1,2}.dirtyDose;
end

if isfield(resultGUI_cell{1,1},'effect') && isfield(resultGUI_cell{1,2},'effect')
    resultGUI.effect = resultGUI_cell{1,1}.effect + resultGUI_cell{1,2}.effect;
elseif isfield(resultGUI_cell{1,1},'effect') && ~isfield(resultGUI_cell{1,2},'effect')
    resultGUI.effect1 = resultGUI_cell{1,1}.effect;
elseif ~isfield(resultGUI_cell{1,1},'effect') && isfield(resultGUI_cell{1,2},'effect')
    resultGUI.effect2 = resultGUI_cell{1,2}.effect;
end

if isfield(resultGUI_cell{1,1},'physicalDose') && isfield(resultGUI_cell{1,2},'physicalDose')
    resultGUI.physicalDose = resultGUI_cell{1,1}.physicalDose + resultGUI_cell{1,2}.physicalDose;
elseif isfield(resultGUI_cell{1,1},'physicalDose') && ~isfield(resultGUI_cell{1,2},'physicalDose')
    resultGUI.physicalDose1 = resultGUI_cell{1,1}.physicalDose;
elseif ~isfield(resultGUI_cell{1,1},'physicalDose') && isfield(resultGUI_cell{1,2},'physicalDose')
    resultGUI.physicalDose2 = resultGUI_cell{1,2}.physicalDose;
end


end