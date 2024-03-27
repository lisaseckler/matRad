function resultGUI = matRad_calcResultGUIstruct(resultGUI,pln)
% creates resultGUI struct by adding the fields in the resultGUI cell

resultGUI_cell = resultGUI;
resultGUI = struct();

if isfield(resultGUI_cell{1,1},'LET') && isfield(resultGUI_cell{1,2},'LET')
    resultGUI.LET = resultGUI_cell{1,1}.LET * pln(1).numOfFractions + resultGUI_cell{1,2}.LET * pln(2).numOfFractions;
% elseif isfield(resultGUI_cell{1,1},'LET') && ~isfield(resultGUI_cell{1,2},'LET')
%     resultGUI.LET1 = resultGUI_cell{1,1}.LET * pln(1).numOfFractions + (ones(size(resultGUI_cell{1,2}.physicalDose)) * 0.3) * pln(2).numOfFractions;
% elseif ~isfield(resultGUI_cell{1,1},'LET') && isfield(resultGUI_cell{1,2},'LET')
%     resultGUI.LET2 = (ones(size(resultGUI_cell{1,1}.physicalDose)) * 0.3) * pln(1).numOfFractions + resultGUI_cell{1,2}.LET * pln(2).numOfFractions;
end

% if isfield(resultGUI_cell{1,1},'RBE') && isfield(resultGUI_cell{1,2},'RBE')
%     resultGUI.RBE = resultGUI_cell{1,1}.RBE * pln(1).numOfFractions + resultGUI_cell{1,2}.RBE * pln(2).numOfFractions;
% elseif isfield(resultGUI_cell{1,1},'RBE') && ~isfield(resultGUI_cell{1,2},'RBE')
%     resultGUI.RBE1 = resultGUI_cell{1,1}.RBE * pln(1).numOfFractions + 1.1 * pln(2).numOfFractions;
% elseif ~isfield(resultGUI_cell{1,1},'RBE') && isfield(resultGUI_cell{1,2},'RBE')
%     resultGUI.RBE2 = 1.1 * pln(1).numOfFractions + resultGUI_cell{1,2}.RBE * pln(2).numOfFractions;
% elseif ~isfield(resultGUI_cell{1,1},'RBE') && ~isfield(resultGUI_cell{1,2},'RBE')
%     resultGUI.RBE3 = 1.1 * pln(1).numOfFractions + 1.1 * pln(2).numOfFractions;
% end

if isfield(resultGUI_cell{1,1},'RBExD') && isfield(resultGUI_cell{1,2},'RBExD')
    resultGUI.RBExD = resultGUI_cell{1,1}.RBExD * pln(1).numOfFractions + resultGUI_cell{1,2}.RBExD * pln(2).numOfFractions;
% elseif isfield(resultGUI_cell{1,1},'RBExD') && ~isfield(resultGUI_cell{1,2},'RBExD')
%     resultGUI.RBExD1 = resultGUI_cell{1,1}.RBExD * pln(1).numOfFractions + resultGUI_cell{1,2}.physicalDose * 1.1 * pln(2).numOfFractions;
% elseif ~isfield(resultGUI_cell{1,1},'RBExD') && isfield(resultGUI_cell{1,2},'RBExD')
%     resultGUI.RBExD2 = resultGUI_cell{1,1}.physicalDose * 1.1 * pln(1).numOfFractions + resultGUI_cell{1,2}.RBExD;
% elseif ~isfield(resultGUI_cell{1,1},'RBExD') && ~isfield(resultGUI_cell{1,2},'RBExD')
%     resultGUI.RBExD3 = resultGUI_cell{1,1}.physicalDose * 1.1 * pln(1).numOfFractions + resultGUI_cell{1,2}.physicalDose * 1.1 * pln(2).numOfFractions;
end

% if isfield(resultGUI_cell{1,1},'SqrtBetaDoseCube') && isfield(resultGUI_cell{1,2},'SqrtBetaDoseCube')
%     resultGUI.SqrtBetaDoseCube = resultGUI_cell{1,1}.SqrtBetaDoseCube + resultGUI_cell{1,2}.SqrtBetaDoseCube;
% elseif isfield(resultGUI_cell{1,1},'SqrtBetaDoseCube') && ~isfield(resultGUI_cell{1,2},'SqrtBetaDoseCube')
%     resultGUI.SqrtBetaDoseCube1 = resultGUI_cell{1,1}.SqrtBetaDoseCube;
% elseif ~isfield(resultGUI_cell{1,1},'SqrtBetaDoseCube') && isfield(resultGUI_cell{1,2},'SqrtBetaDoseCube')
%     resultGUI.SqrtBetaDoseCube2 = resultGUI_cell{1,2}.SqrtBetaDoseCube;
% end
% 
% if isfield(resultGUI_cell{1,1},'alpha') && isfield(resultGUI_cell{1,2},'alpha')
%     resultGUI.alpha = resultGUI_cell{1,1}.alpha + resultGUI_cell{1,2}.alpha;
% elseif isfield(resultGUI_cell{1,1},'alpha') && ~isfield(resultGUI_cell{1,2},'alpha')
%     resultGUI.alpha1 = resultGUI_cell{1,1}.alpha;
% elseif ~isfield(resultGUI_cell{1,1},'alpha') && isfield(resultGUI_cell{1,2},'alpha')
%     resultGUI.alpha2 = resultGUI_cell{1,2}.alpha;
% end
% 
% if isfield(resultGUI_cell{1,1},'alphaDoseCube') && isfield(resultGUI_cell{1,2},'alphaDoseCube')
%     resultGUI.alphaDoseCube = resultGUI_cell{1,1}.alphaDoseCube + resultGUI_cell{1,2}.alphaDoseCube;
% elseif isfield(resultGUI_cell{1,1},'alphaDoseCube') && ~isfield(resultGUI_cell{1,2},'alphaDoseCube')
%     resultGUI.alphaDoseCube1 = resultGUI_cell{1,1}.alphaDoseCube;
% elseif ~isfield(resultGUI_cell{1,1},'alphaDoseCube') && isfield(resultGUI_cell{1,2},'alphaDoseCube')
%     resultGUI.alphaDoseCube2 = resultGUI_cell{1,2}.alphaDoseCube;
% end
% 
% if isfield(resultGUI_cell{1,1},'beta') && isfield(resultGUI_cell{1,2},'beta')
%     resultGUI.beta = resultGUI_cell{1,1}.beta + resultGUI_cell{1,2}.beta;
% elseif isfield(resultGUI_cell{1,1},'beta') && ~isfield(resultGUI_cell{1,2},'beta')
%     resultGUI.beta1 = resultGUI_cell{1,1}.beta;
% elseif ~isfield(resultGUI_cell{1,1},'beta') && isfield(resultGUI_cell{1,2},'beta')
%     resultGUI.beta2 = resultGUI_cell{1,2}.beta;
% end
% 
% if isfield(resultGUI_cell{1,1},'cleanDose') && isfield(resultGUI_cell{1,2},'cleanDose')
%     resultGUI.cleanDose = resultGUI_cell{1,1}.cleanDose + resultGUI_cell{1,2}.cleanDose;
% elseif isfield(resultGUI_cell{1,1},'cleanDose') && ~isfield(resultGUI_cell{1,2},'cleanDose')
%     resultGUI.cleanDose1 = resultGUI_cell{1,1}.cleanDose;
% elseif ~isfield(resultGUI_cell{1,1},'cleanDose') && isfield(resultGUI_cell{1,2},'cleanDose')
%     resultGUI.cleanDose2 = resultGUI_cell{1,2}.cleanDose;
% end

if isfield(resultGUI_cell{1,1},'dirtyDose') && isfield(resultGUI_cell{1,2},'dirtyDose')
    resultGUI.dirtyDose = resultGUI_cell{1,1}.dirtyDose * pln(1).numOfFractions + resultGUI_cell{1,2}.dirtyDose * pln(2).numOfFractions;
% elseif isfield(resultGUI_cell{1,1},'dirtyDose') && ~isfield(resultGUI_cell{1,2},'dirtyDose')
%     resultGUI.dirtyDose1 = resultGUI_cell{1,1}.dirtyDose * pln(1).numOfFractions;
% elseif ~isfield(resultGUI_cell{1,1},'dirtyDose') && isfield(resultGUI_cell{1,2},'dirtyDose')
%     resultGUI.dirtyDose2 = resultGUI_cell{1,2}.dirtyDose * pln(2).numOfFractions;
end

% if isfield(resultGUI_cell{1,1},'effect') && isfield(resultGUI_cell{1,2},'effect')
%     resultGUI.effect = resultGUI_cell{1,1}.effect + resultGUI_cell{1,2}.effect;
% elseif isfield(resultGUI_cell{1,1},'effect') && ~isfield(resultGUI_cell{1,2},'effect')
%     resultGUI.effect1 = resultGUI_cell{1,1}.effect;
% elseif ~isfield(resultGUI_cell{1,1},'effect') && isfield(resultGUI_cell{1,2},'effect')
%     resultGUI.effect2 = resultGUI_cell{1,2}.effect;
% end

if isfield(resultGUI_cell{1,1},'physicalDose') && isfield(resultGUI_cell{1,2},'physicalDose')
    resultGUI.physicalDose = resultGUI_cell{1,1}.physicalDose * pln(1).numOfFractions + resultGUI_cell{1,2}.physicalDose * pln(2).numOfFractions;
% elseif isfield(resultGUI_cell{1,1},'physicalDose') && ~isfield(resultGUI_cell{1,2},'physicalDose')
%     resultGUI.physicalDose1 = resultGUI_cell{1,1}.physicalDose * pln(1).numOfFractions;
% elseif ~isfield(resultGUI_cell{1,1},'physicalDose') && isfield(resultGUI_cell{1,2},'physicalDose')
%     resultGUI.physicalDose2 = resultGUI_cell{1,2}.physicalDose * pln(2).numOfFractions;
end

% if isfield(resultGUI_cell{1,1},'mLETDose') && isfield(resultGUI_cell{1,2},'mLETDose')
%     resultGUI.physicalDose = resultGUI_cell{1,1}.physicalDose + resultGUI_cell{1,2}.physicalDose;
% elseif isfield(resultGUI_cell{1,1},'mLETDose') && ~isfield(resultGUI_cell{1,2},'mLETDose')
%     resultGUI.physicalDose1 = resultGUI_cell{1,1}.physicalDose;
% elseif ~isfield(resultGUI_cell{1,1},'mLETDose') && isfield(resultGUI_cell{1,2},'mLETDose')
%     resultGUI.physicalDose2 = resultGUI_cell{1,2}.physicalDose;
% end


end