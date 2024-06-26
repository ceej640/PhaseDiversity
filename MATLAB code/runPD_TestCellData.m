function runPD_TestCellData(Cycle)
%Run script to replicate on-line call of AO data for U2OS cell data
%INPUT: Cycle number. %1 = Measure Aberrated WF, 2...10 = Corrections n-1.
pkgdir = findpkg();
AOstruct = getTestInputs(768);
AOpath = fullfile(pkgdir, '/Datasets/AO/230921 AO0057 U2OS_Cell');

[coeffsOut, ~, runtime, imgEstimate, waveFront, wRMSE, imageQualityMetric] = processPhaseDiversityImages(AOpath,['Iteration ' num2str(Cycle)], AOstruct.fname, ...
    AOstruct.nImages, AOstruct.nRepetitions, AOstruct.nZernike, AOstruct.iterationLimit, AOstruct.gamma, 0, AOstruct.cropSize, AOstruct.backgroundValue,...
    AOstruct.flagShowInput, AOstruct.flagShowRecon, detectGPU(AOstruct.flagGPU), 'ANSI', 1, 0, 1);

function matchedDir = findpkg()
    pathdir = path;
    targetDir = [filesep 'PhaseDiversity'];
    if contains(pathdir, targetDir)
        pathDirs = strsplit(path, ';');
        matchedDir = '.';
        for i = 1:length(pathDirs)
            if endsWith(pathDirs{i}, targetDir)
                matchedDir = pathDirs{i};
                break;
            end
        end
    end
end

function AOstruct = getTestInputs(Sx)
    AOstruct = struct('fname','Stack','nImages',5,'nRepetitions',1,'nZernike',20,'iterationLimit',100,'gamma',1e-6,'cropSize',Sx,'backgroundValue',400,'flagShowInput',0, ...
        'flagShowRecon',0,'flagGPU','auto');
end
end