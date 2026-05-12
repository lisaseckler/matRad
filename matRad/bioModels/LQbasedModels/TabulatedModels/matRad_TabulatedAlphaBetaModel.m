classdef matRad_TabulatedAlphaBetaModel < matRad_TabulatedDoseAveragedKernelModel

    properties (Constant)
        model      = 'doseAveragedTabulatedAlphaBeta'
        quantitiesToAverage = {'alpha', 'sqrtBeta', 'zs'};

        quantitiesInTable = {'alpha', 'beta', 'zs'};
        requiredQuantities = {'Fluence'};
        possibleRadiationModes = {'protons', 'carbon', 'helium'};

        cellType = ["HSG","V79","CHO"];
    
    end

    methods
        function this = matRad_TabulatedAlphaBetaModel()
            this@matRad_TabulatedDoseAveragedKernelModel;
        end

        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel,kernels,cstDoseGrid)
            % This function assignis the interpolated kernels to the
            % correct tissue class
            
            bixel  = calcBiologicalQuantitiesForBixel@matRad_LQBasedModel(this,bixel);
            bixel.sqrtBeta = NaN*ones(numel(bixel.radDepths),1);
            % Assume the kernels for multiple quantities for the same model
            % have the same dimension
            if ~isempty(this.quantityTable) && isempty(this.ZstarTable)
                if ~isfield(this.quantityTable,'zs')
                    correctedQuantities = this.quantitiesToAverage;
                    correctedQuantities(strcmp(correctedQuantities,'zs')) = [];
                end
            else
                if ~isfield(this.ZstarTable,'alpha') && ~isfield(this.ZstarTable,'beta')
                    correctedQuantities = this.quantitiesToAverage;
                    correctedQuantities(strcmp(correctedQuantities,'alpha')) = [];
                    correctedQuantities(strcmp(correctedQuantities,'sqrtBeta')) = [];
                end
            end

            numOfTissueClass = size(bixel.baseData.(correctedQuantities{1}),2);
            
            for i = 1:numOfTissueClass
                mask = bixel.vTissueIndex == i;
  
                if any(mask)
                    
                    for quantityName = correctedQuantities
                        bixel.(quantityName{1})(mask) = kernels.(quantityName{1})(mask,i);
                        if ismember('zs',quantityName)
                            bixel.beta(:) = cstDoseGrid{1,5}.bioParams.betaX;
                        end
                    end
                end
            end       
        end

        function outQuantity = interpolateQuantityOnSpectra(this, spectra)
            outQuantity = interpolateQuantityOnSpectra@matRad_TabulatedDoseAveragedKernelModel(this,spectra);
            if ~isempty(this.quantityTable)
                for i=1:numel(outQuantity)
                    outQuantity(i).beta(outQuantity(i).beta<0) = 0; % so that we don't get a complex sqrtBeta with the MCF machine
                    outQuantity(i).sqrtBeta = sqrt(outQuantity(i).beta);
                end
            end
        end

        function doseAveragedKernels = computeDoseAveragedKernels(this,quantity, spectra, sp)
            doseAveragedKernels = computeDoseAveragedKernels@matRad_TabulatedDoseAveragedKernelModel(this,quantity, spectra, sp);
            if ~isempty(this.quantityTable)
                if isfield(this.quantityTable.meta.modelParameters,'cellLine')
                    doseAveragedKernels.cellLine = this.quantityTable.meta.modelParameters.cellLine;
                elseif isfield(this.quantityTable.data,'cellLine')
                    doseAveragedKernels.cellLine = this.quantityTable.data(1).cellLine;
                end
            elseif ~isempty(this.ZstarTable)
                 if isfield(this.ZstarTable.meta.modelParameters,'cellLine')
                    doseAveragedKernels.cellLine = this.ZstarTable.meta.modelParameters.cellLine;
                 end
            end
        end

        function tableData = getTableDataForAlphaBeta(this, alphaX, betaX,cellLine)
            % Scroll the RBE table and check if the requested alpha beta
            % ratio is available

            tableData = [];
            if isempty(this.quantityTable)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('RBEtable not found');
            elseif isfield(this.quantityTable.meta,'alphaX') && isfield(this.quantityTable.meta,'betaX')
                for i=1:numel(this.quantityTable.meta.alphaX)
                    if (this.quantityTable.meta.alphaX(i) == alphaX) && (this.quantityTable.meta.betaX(i) == betaX)
                        tableData = this.extractFragmentsWithZA([this.includedFragments.Z], [this.includedFragments.A], this.quantityTable.data);
                         
                        for j=1:numel(tableData)
                            tableData(j).alpha = tableData(j).alpha(:,i);
                            tableData(j).beta  = tableData(j).beta(:,i);
                        end

                        continue;
                    end
                end
                if isempty(tableData)
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning(sprintf('No match found for alpha/beta = %f/%f in RBE table', alphaX, betaX));                    
                end
            elseif isfield(this.quantityTable.meta.modelParameters,'cellLine')
                for i=1:numel(this.quantityTable.meta.modelParameters.cellLine)
                    if (this.quantityTable.meta.modelParameters.cellLine(i) == cellLine)
                        tableData = this.extractFragmentsWithZA([this.includedFragments.Z], [this.includedFragments.A], this.quantityTable.data);
                         
                        for j=1:numel(tableData)
                            tableData(j).cellLine = cellLine(:,i);
                        end

                        continue;
                    end
                end
            end
        end

    end

    methods (Static)
        function RBEtable = loadTable(fileName)

            % This function loads the specified RBE table
            matRad_cfg = MatRad_Config.instance();

            searchPath = {fullfile(matRad_cfg.matRadSrcRoot,'bioModels','RBEtables'),...    % default matrad folder
                fullfile(matRad_cfg.primaryUserFolder, 'RBEtables')};             % user defined RBE table

            try
                load(fullfile(searchPath{1}, [fileName, '.mat']), 'RBEtable');

            catch
                try
                    load(fullfile(searchPath{2}, [fileName, '.mat']), 'RBEtable');

                catch
                    matRad_cfg.dispError('Cannot find RBEtable: %s', fileName);
                end
            end
        end

    end
end