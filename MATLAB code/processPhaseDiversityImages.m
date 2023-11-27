function [coeffsOut, coeffRMSE, runtime, imgEstimate, waveFront, wRMSE, imageQualityMetric] = processPhaseDiversityImages(path, iterationFolder, filename, ...
    nImages, nRepetitions, nZernike, iterationLimit, gamma, alpha, cropSize, backgroundValue,...
    flagShowInput, flagShowRecon, flagGPU, zernikeIndexing, penalChoice, GraduatedOptimization, BinFactor)
% By Courtney Johnson, Magdalena Schneider and Min Guo
%16OCT23
%Use with Dataset: XXXXXXXXXXXXXXXXXXXXX
%Wrapper function for on-line wavefront estimation: Script takes diversity image stacks and associated Zernike coeffs., and estimates wavefront and object.
% 
% 
% OUTPUT:
% 
% coeffsOut: reconstructed Zernike coefficients, to be applied to DM to cancel aberration
% 
% coeffRMSE: root mean squared error of coefficient estimate
% 
% runtime: time to estimate wavefront
% 
% imgEstimate: Image output of Estimated Object
% 
% waveFront: Estimated wavefront
% 
% wRMSE: root mean squared error of wavefront estimate
% 
% imageQualityMetric: DCTS, a measure of image sharpness
%
% 
% INPUT:
% 
% path: character vector of path of folder containing all cycles
% 
% iterationFolder: Subfolder character vector. Location of Image files. 
% ***Here iteration refers to the cycle number, where 1 is the aberrated, uncorrected image. ie: 'Iteration 1'***
% 
% filename: file name of input image stack without extension as character vector. Typically 'Stack'
% 
% nImages: total number of images (raw image + diversity images)
% 
% nRepetitions: number of images to be averaged (Typically 1)
% 
% nZernike: Maximum Zernike coefficient to be estimated
% 
% iterationLimit: maximum number of loops for the reconstruction algorithm
% 
% gamma: tuning parameter for image penalty term, typically set at 1e-6
% 
% alpha: set to 0, not used.
% 
% cropSize: image size to crop to for calculation. Crop location is centered within
% image frame.
% 
% backgroundValue: value to subtract from all pixels
% 
% flagShowInput: 1 = output figures showing input diversity phases
% 
% flagShowRecon: 1 = output figures showing result
% 
% flagGPU: 1 = use GPU, 0 = use CPU
% 
% zernikeIndexing: Character vector specifying ordering of Zernike coeffs. Typically 'ANSI'
% 
% penalChoice: Specify penalty term in algorithm,  Typically not used, set to 1 = penalize square norm
% 
% GraduatedOptimization: Estimate only lower order Zernike coeffs first. Set to 0, not used.
% 
% BinFactor: Bin image pixels and scale pixel size accordingly to reduce
% image size. Typically not used, set to 1.
% 

plotFigures=false;
tempname = strsplit(strip(fullfile(path, iterationFolder),filesep));
ItNumber = str2double(tempname{end});
tic

%% Setup parameters

if nargin == 15
    penalChoice = 1; % default for Tikhonov norm
end

pixelSize = 0.104*BinFactor; % um
wavelength = 0.532; % um
NA = 1.2;
rotAng = 0;

imageFilepath = fullfile(path, iterationFolder, [filename '.tif']); %Location of Image stack
zernikeCoefficientsFilepath = fullfile(path,'CoeffsIn.txt'); %Location of text file containing known Zernike coeffs.

%% Create output folders
iterationDir = fullfile(path, iterationFolder);
if isequal(exist(iterationDir, 'dir'),7)
    disp(['output folder:' iterationDir]);
else
    mkdir(iterationDir);
    disp(['output folder created:' iterationDir]);
end

savedir = fullfile(path, 'results');
if ~exist(savedir, 'dir')
    mkdir(savedir)
end


switch zernikeIndexing
    case 'ANSI'
        conType = 'ANSI2ANSI';
    otherwise
        error('Wrong index type! Supports ANSI only.');
end

% Phases and zernike coefficients
% Exclude from estimation 0: piston; 1: tilt Y; 2: tilt X;
zernikeIndices = 3:nZernike;


%% Preprocess images

% Read images
disp('Preprocessing images...');
images = fileIO_lvtiff2mat(imageFilepath, nImages, nRepetitions, cropSize, backgroundValue);
WriteTifStack(images, fullfile(path, iterationFolder, 'Image_phasediversity.tif'), 32);

% Set Zernike coefficients to be estimated
if GraduatedOptimization==1
    % upboosting steps, e.g., 4th (14), 5th(20)
    oIdx = [10 14 20];
    if nZernike<=oIdx(1)
        zernSteps = nZernike - zernikeIndices(1) + 1;
    elseif nZernike<=oIdx(2)
        zernSteps = [oIdx(1) nZernike] - zernikeIndices(1) + 1;
    elseif nZernike<=oIdx(3)
        zernSteps = [oIdx(1) oIdx(2) nZernike] - zernikeIndices(1) + 1;
    else
        error('processPhaseDiversityImages: Zernike coefficents order out of range');
    end
else
    zernSteps = nZernike - zernikeIndices(1) + 1;
end


% Read Zernike coefficients from file
coeffsRaw = importdata(zernikeCoefficientsFilepath);
if size(coeffsRaw, 2) > nZernike
    coeffsRaw = coeffsRaw(:,1:nZernike);
end
if size(coeffsRaw, 1) < nImages % number of applied phase diversities
    error('Phase number is less than image number');
end

% Set range of Zernike indices to fit
minZernikeIndex = 3; % 1: tilt Y; 2: tilt X; 3: defocus;
maxZernikeIndex = nZernike;

% Unknown aberration
coeffsInitial = coeffsRaw(1, minZernikeIndex:maxZernikeIndex);
if ItNumber~=1 % Target Coeffs for correction are 0, not test aberration.
    coeffsInitial=zeros(1,length(zernikeIndices));
end

% Phase diversity
coeffs_div = coeffsRaw(2:nImages, minZernikeIndex:maxZernikeIndex);
coeffs_all = [coeffsInitial; coeffs_div];


%% Display runtime
runtime1 = toc;

%% Reconstruct Zernike coefficients

[cEstimate, imgEstimate, ~] = reconstructZernikeAberrations(images, zernikeIndices, coeffs_div, gamma,...
    alpha, iterationLimit, zernSteps, pixelSize, wavelength, NA, flagGPU, penalChoice, rotAng, savedir,plotFigures);

runtime2 = toc;

imgEstimate = double(rescale(imgEstimate));
coeffRMSE = double(rms(cEstimate(:) - coeffsInitial(:)));

% Display runtime

disp(['... Time cost: ', num2str(runtime2-runtime1)]);


% Save



%% Calculate PSFs and waveFront based on input zernike coefficients

baseImage = squeeze(images(:,:,1));
[~, PSFs, waveFronts] = createPhaseDiversityImages(baseImage, zernikeIndices, coeffs_all, pixelSize, wavelength, NA, 'none', 10, rotAng, flagGPU);
WriteTifStack(PSFs, fullfile(path, iterationFolder, 'PSF_phasediversity.tif'), 32);
WriteTifStack(imgEstimate, fullfile(path, iterationFolder, 'Image_estimated.tif'), 32);

% Estimated wavefront
[Sx, Sy] = size(baseImage);
parMask.nGrid = Sx;
parMask.rotationAngle = deg2rad(rotAng);
parMask.maxRadius = calculatePupilRadius(Sx, pixelSize, wavelength, NA);
pupilMask = Mask(parMask);
idx = logical(pupilMask.values);
zernikeObject = ZernikePolynomials.getInstance(pupilMask,'ANSI',0,flagGPU);
waveFront = double(gather(zernikeObject.getAberration(zernikeIndices, cEstimate) .* idx));
wRMSE = gather(rms(rms(waveFront(idx))));


%% Assess image quality

imageQualityMetric = calculateImageQualityMetrics(baseImage, 1);

%% Export zernike coefficient to txt file: imagine optic convention

fileTxtOutname = fullfile(path, iterationFolder, 'zernCoeffs_estimated.txt');
coeffsOut=zeros(1,maxZernikeIndex);
coeffsOut(zernikeIndices)=cEstimate;
dlmwrite(fileTxtOutname, coeffsOut ,'delimiter','\t','precision','%.6f');


%% Plot results

if flagShowInput == 1
    a = 20; % size of PSF images for show
    F1 = figure; % input images and phases

    figure(F1); subplot(nImages,3,1);
    plotWaveFront(waveFronts(:,:,1))
    title('Wavefront: aberrated');

    figure(F1); subplot(nImages,3,2);
    PSF = PSFs(:,:,1);
    PSF = PSF/max(PSF(:));
    PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
    plotImage(PSF)
    title('PSF:aberrated');

    figure(F1); subplot(nImages,3,3);
    img = images(:,:,1);
    img = img/max(img(:));
    plotImage(img)
    title('Image:aberrated');

    figure(F1); subplot(nImages,3,4);
    plotWaveFront(waveFronts(:,:,2))
    title('Add:phase1');

    figure(F1); subplot(nImages,3,5);
    PSF = PSFs(:,:,2);
    PSF = PSF/max(PSF(:));
    PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
    plotImage(PSF)
    title('PSF:phase1');

    figure(F1); subplot(nImages,3,6);
    plotImage(images(:,:,2)/max(img(:)))
    title('Image:phase1');

    if nImages >= 3
        for i = 3:nImages
            figure(F1); subplot(nImages,3,3*(i-1)+1);
            plotWaveFront(waveFronts(:,:,i))
            title(['Add:phase', num2str(i-1)]);

            figure(F1); subplot(nImages,3,3*(i-1)+2);
            PSF = PSFs(:,:,i);
            PSF = PSF/max(PSF(:));
            PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
            plotImage(PSF)
            title(['PSF:phase', num2str(i-1)]);

            figure(F1); subplot(nImages,3,3*(i-1)+3);
            img = images(:,:,i);
            img = img/max(img(:));
            plotImage(img)
            title(['Image:phase', num2str(i-1)]);
        end
    end
    savefig(fullfile(path, iterationFolder, 'input.fig'));
end


if flagShowRecon == 1
    [~, ~, idx] = def_pupilcoor(Sx, pixelSize, wavelength, NA);
    waveFront_original = waveFronts(:,:,1);

    wMin = min(waveFront_original(:));
    wMax = max(waveFront_original(:));
    if wMin == wMax
        wMin = -5;
        wMax = 5;
    end

    waveFrontForShow = nan(Sx,Sx);
    waveFrontForShow(idx) = waveFront_original(idx);

    F2 = figure; subplot(2,2,1);
    plotWaveFront(waveFrontForShow)
    clim([wMin wMax])
    title('Wavefront: aberrated');

    figure(F2); subplot(2,2,2);
    plotImage(baseImage/max(baseImage(:)))
    title('Image: aberrated');

    figure(F2); subplot(2,2,3);
    waveFrontForShow = nan(Sx,Sx);
    waveFrontForShow(idx) = waveFront(idx);
    plotWaveFront(waveFrontForShow)
    clim([wMin wMax])
    title('Wavefront: estimated');

    figure(F2); subplot(2,2,4);
    plotImage(imgEstimate/max(imgEstimate(:)))
    title('Image:estimated');
    savefig(fullfile(path, iterationFolder, 'retrieval.fig'));
    
    % Plot Zernike coefficients as MATLAB convention
    figure
    plotZernikeCoefficients(zernikeIndices, coeffsInitial, cEstimate)
    title('Zernike Modes (MATLAB)');
    savefig(fullfile(path, iterationFolder, 'coeff.fig'));
    
end

runtime=runtime2-runtime1;
save(fullfile(path, iterationFolder, 'data.mat'));

end


% Functions for plotting

function plotWaveFront(wavefront)
    xi = 1:size(wavefront, 1);
    pcolor(xi,xi,wavefront), shading interp
    axis square, Fc = colorbar;
    xlabel(Fc,'\mum');
end

function plotImage(image)
    imshow(image,[]),colorbar;
end

function plotZernikeCoefficients(zernikeIndices, coeff, estimate)
    plot(zernikeIndices, coeff, zernikeIndices, estimate,'LineWidth',2);
    legend( 'target','estimated');
    xlabel('Zernike Index');
    ylabel('Coeff Value');
    set(gca,'FontSize', 14);
end