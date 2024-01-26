classdef matRad_MaxDVHDirtyDose < DirtyDoseObjectives.matRad_DirtyDoseObjective
% matRad_MaxDVHDirtyDose Implements a penalized maximum dirtyDose DVH objective
%   See matRad_DirtyDoseObjective for interface description
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
        name = 'Max DVH';
        parameterNames = {'d', 'V^{max}'};
        parameterTypes = {'dirtyDose','numeric'};
    end
    
    properties
        parameters = {30,95};
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MaxDVHDirtyDose(penalty,dRef,vMaxPercent)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DirtyDoseObjectives.matRad_DirtyDoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin >= 3 && isscalar(vMaxPercent)
                    obj.parameters{2} = vMaxPercent;
                end
                
                if nargin >= 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end        
        
        %% Calculates the Objective Function value
        function fDirtyDose = computeDirtyDoseObjectiveFunction(obj,dirtyDose)                       
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = dirtyDose - obj.parameters{1};

            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVHDirtyDose(refVol,dirtyDose);

            
            deviation(dirtyDose < obj.parameters{1} | dirtyDose > d_ref2) = 0;
   
            % calculate objective function
            fDirtyDose = (1/numel(dirtyDose))*(deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fDirtyDoseGrad   = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = dirtyDose - obj.parameters{1};
            
            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVHDirtyDose(refVol,dirtyDose);
            
            deviation(dirtyDose < obj.parameters{1} | dirtyDose > d_ref2) = 0;

            % calculate delta
            fDirtyDoseGrad = (2/numel(dirtyDose))*deviation;
        end
    end
    
end