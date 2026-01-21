function Data =  loadSimuDataDose(foldername, resolution,E, dataType)


    foldername = [foldername, sprintf('Energy%g',E), '\Results\']; 
    switch dataType
        case 'binary'
            filename = sprintf('data_ED_Coarse.bin');
            filename = [foldername, filename];
            fID = fopen(filename);
            DataTmp = fread(fID,resolution.zCoarseNum*resolution.rNum,'double');
            fclose(fID);
            DataTmp = reshape(DataTmp, [resolution.rNum,resolution.zCoarseNum,]); 
            Data.ED_Coarse = DataTmp'./resolution.numPrimary./resolution.zCoarse;

            filename = sprintf('data_ED_Fine.bin');
            filename = [foldername, filename];
            fID = fopen(filename);
            DataTmp = fread(fID,resolution.zFineNum*resolution.rNum,'double');
            fclose(fID);
            DataTmp = reshape(DataTmp, [resolution.rNum,resolution.zFineNum,]); 
            Data.ED_Fine = DataTmp'./resolution.numPrimary./resolution.zFine;
        otherwise

            error('DataType not implemented')
    end
    fclose all;


end

    

    
