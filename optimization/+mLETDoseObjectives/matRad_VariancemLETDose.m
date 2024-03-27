classdef matRad_VariancemLETDose < mLETDoseObjectives.matRad_mLETDoseObjective
% matRad_VariancemLETDose Implements a homogenous dose distribution
%   See matRad_mLETDoseObjective for interface description
%
% References 
%     -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
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
        name = 'Variance';
        parameterNames = {'dRef'};
        parameterTypes = {'mLETDose'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end
    
    methods
        function obj = matRad_VariancemLETDose(penalty)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@mLETDoseObjectives.matRad_mLETDoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
   
        %% Calculates the Objective Function value
        function fmLETDose = computemLETDoseObjectiveFunction(obj,mLETDose)
            % deviation : dose minus prefered dose
            %deviation = dose - obj.parameters{1};
            %fClusterDose = clusterDose'*clusterDose / numel(clusterDose) - mean(clusterDose)^2;
            %fClusterDose = fClusterDose * obj.penalty*numel(clusterDose)/(numel(clusterDose) - 1);
            fmLETDose = obj.penalty * var(mLETDose);
            
            % claculate objective function
            %fClusterDose = obj.penalty/numel(clusterDose) * (deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fmLETDoseGrad   = computemLETDoseObjectiveGradient(obj,mLETDose)
            % deviation : Dose minus prefered dose
            % deviation = dose - obj.parameters{1};
            %deviation = (clusterDose - obj.parameters{1});
            
            % calculate delta
            %fClusterDoseGrad = 2 * obj.penalty/numel(clusterDose) * deviation;
            fmLETDoseGrad = obj.penalty * 2/(numel(mLETDose) - 1) * (mLETDose - mean(mLETDose));
        end
    end
    
end