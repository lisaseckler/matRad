classdef matRad_SquaredUnderdosingLETt < LETtObjectives.matRad_LETtObjective
% matRad_SquaredUnderdosingLETt Implements a penalized squared underdosing LETt objective
%   See matRad_LETtObjectives for interface description
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        name = 'Squared Underdosing LETt';
        parameterNames = {'LETt^{min}'};
        parameterTypes = {'LETt'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredUnderdosingLETt(penalty,LETtMin)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@LETtObjectives.matRad_LETtObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin >= 2 && isscalar(LETtMin)
                    obj.parameters{1} = LETtMin;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fLETt = computeLETtObjectiveFunction(obj,LETt)
            % underLETt : LETt minus prefered LETt
            underLETt = LETt - obj.parameters{1};
            
            % apply positive operator
            underLETt(underLETt>0) = 0;
            
            % calculate objective function
            fLETt = 1/numel(LETt) * (underLETt'*underLETt);
        end
        
        %% Calculates the Objective Function gradient
        function fLETtGrad   = computeLETtObjectiveGradient(obj,LETt)
            % underLETt : LETt minus prefered LETt
            underLETt = LETt - obj.parameters{1};
            
            % apply positive operator
            underLETt(underLETt>0) = 0;
            
            % calculate delta
            fLETtGrad = 2/numel(LETt) * underLETt;
        end
    end
    
end