function gpuFlag = detectGPU(gpuFlag)
    arguments
        gpuFlag {mustBeMember(gpuFlag,[0,1,'auto'])}   
    end
    if strcmp(gpuFlag, 'auto')
        if gpuDeviceCount("available")>0
            gpuFlag = 1;
            fprintf('Running on GPU...\n')
        else
            gpuFlag = 0;
            fprintf('Running on CPU...\n')
        end
    end
end