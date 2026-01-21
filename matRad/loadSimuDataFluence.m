function Data =  loadSimuDataFluence(foldername, params,E, dataType)

    foldername = [foldername, sprintf('Energy%g',E), '\Results\Fluence\']; 
    Eres = params.spectrum_ebinsMax/params.spectrum_ebinsNum;
    Data.energyBin = linspace(Eres/2, params.spectrum_ebinsMax -Eres/2, params.spectrum_ebinsNum);
    Eres_electron = params.electron_spectrum_ebinsMax/params.electron_spectrum_ebinsNum;
    Data.energyBinEl = linspace(Eres_electron /2, params.electron_spectrum_ebinsMax - Eres_electron/2, params.electron_spectrum_ebinsNum);
    Data.Z = [1,1,1,2,3,4,5,6,7,8,-1];
    Data.A = [1,2,3,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
    switch dataType
        case 'binary'
            for i = 1:numel(Data.Z)-1
                if isnan(Data.A(i))
                    filename =  [foldername, sprintf('Fluence_Z%i_phantom_1.bin',Data.Z(i))]; 
                else
                    filename =  [foldername, sprintf('Fluence_Z%i_A%i_phantom_1.bin',Data.Z(i),Data.A(i))]; 
                end
                Data.spectrum(i).Coarse = loadData(filename,params.spectrum_ebinsNum,params.zCoarseNum)./params.numPrimary;
                if isnan(Data.A(i))
                    filename =  [foldername, sprintf('Fluence_Z%i_phantom_2.bin',Data.Z(i))]; 
                else
                    filename =  [foldername, sprintf('Fluence_Z%i_A%i_phantom_2.bin',Data.Z(i),Data.A(i))]; 
                end
                Data.spectrum(i).Fine =  loadData(filename,params.spectrum_ebinsNum,params.zFineNum)./params.numPrimary;
            end
            filename =  [foldername, sprintf('Fluence_electrons_phantom_1.bin')]; 
            Data.spectrum(numel(Data.Z)).Coarse = loadData(filename,params.electron_spectrum_ebinsNum,params.zCoarseNum)./params.numPrimary;
            filename =  [foldername, sprintf('Fluence_electrons_phantom_2.bin')];
            Data.spectrum(numel(Data.Z)).Fine = loadData(filename,params.electron_spectrum_ebinsNum,params.zFineNum)./params.numPrimary;

        otherwise

            error('DataType not implemented')
    end
    fclose all;


end

function result = loadData(filename,nEBins,nZBins)
    fID = fopen(filename);
    result = fread(fID,(nEBins+2)*nZBins, 'double');
    result = reshape(result, [nEBins+2,nZBins]);
    result = result(2:end-1,:);
    fclose(fID);

end
    

    
