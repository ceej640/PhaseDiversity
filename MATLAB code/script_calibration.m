%% Phase Diversity-based DM Calibration Processing Code
% By Courtney Johnson, Magdalena Schneider and Min Guo
%16OCT23
%Use with Dataset: XXXXXXXXXXXXXXXXXXXXX
%For DM Calibration: Script takes defocus image stacks obtained from DM calibration procedure, estimates wavefronts, and returns command matrix.

clear
close all

expID = 'Test'; %Subfolder Label for Saving Outputs

SaveFunctionFigures = false; %false when performing speed tests
SaveScriptFigures = false; %false when performing speed tests
SaveMAT = true; %false when performing speed tests
SaveCM = true; %Save command matrix
CMSH = false; % Rotate command matrix output orientation for validations using Shack-Hartmann

%% System Parameters

pixelSize = 0.104; %Image Pixel Size
lambda = 0.532; % Central Emission Wavelength, um
NA = 1.2; %Numerical Aperture
Sx = 128; %Image Crop size to be operated on
Sy = Sx; %Square image
numz = 5; %Number of slices in image stack
zstepsize = 1; %Space between slices in um
bgval = 500; %Subtract constant bg val

%% Reconstruction settings

flagGPU = 'auto'; % 'auto', 1, 0
zernikeIndicesForFitting = 3:20; % Exclude 0: piston; 1:tilt X; 2: tilt Y from algorithm estimation;
zernikeIndicesIncludingPistonAndTilts = 1:20; % polynomials with tip/tilt
zernNum = length(zernikeIndicesForFitting);
itLimit = 100; % max iteration number
gamma = 1e-6;
zernSteps = zernikeIndicesForFitting(zernNum) - zernikeIndicesForFitting(1) + 1;

%% Initiation

invthresh=0.005;
ft = fittype('poly1');
HSx = Sx/2;
HSy = Sy/2;
xi = 1:Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
RunTime = nan(312,1);
[xpup, ypup] = ind2sub(size(idx),find(idx));
minx = min(xpup(:));
maxx = max(xpup(:));
miny = min(ypup(:));
maxy = max(ypup(:));
rotAng = 90; % Orientation difference of Shack-Hartmann
ds = -zstepsize*(numz-1)/2:zstepsize:zstepsize*(numz-1)/2;  % defocus steps in um
ii = 0;

% GPU
flagGPU = detectGPU(flagGPU);

% Pupil and Zernikes
parMask.nGrid = Sx;
parMask.rotationAngle = deg2rad(0); % rotation is applied later
parMask.maxRadius = calculatePupilRadius(Sx, pixelSize, lambda, NA);
pupilMask = Mask(parMask);
idx = logical(pupilMask.values);
zernikeObject = ZernikePolynomials.getInstance(pupilMask,'ANSI',0,flagGPU);


%% Get list of desired folders

maindir =uigetdir('','Select "DM Calibration" folder in PhaseDiversity/Datasets/');
cd(maindir)
dirlistmain = dir(maindir);
% Filter subfolders
dirlistmain = dirlistmain([dirlistmain(:).isdir]);
% Remove '.' and '..'
dirlistmain = dirlistmain(~ismember({dirlistmain(:).name},{'.','..'}));
% Remove any folder including 'bg' in name
dirlistmain = dirlistmain(~contains({dirlistmain(:).name},{'bg'}));
% Remove any folder including 'Compare Calibration Results'
dirlistmain = dirlistmain(~contains({dirlistmain(:).name},{'Compare Calibration Results'}));

%% Loop

for folderIndex = 1:length(dirlistmain)

    tStart = tic;

    %% Prepare folders

    fileFolderIn = [dirlistmain(folderIndex).folder '/' dirlistmain(folderIndex).name];
    fileFolderPSF = [fileFolderIn '/Processed Cal Data_' expID '/PSF_all/'];
    fileFolderWavefront = [fileFolderIn '/Processed Cal Data_' expID '/wavefront_all/'];

    if SaveCM || SaveFunctionFigures || SaveScriptFigures
        if ~exist(fileFolderPSF, 'dir')
            mkdir(fileFolderPSF);
        end
        if ~exist(fileFolderWavefront, 'dir')
            mkdir(fileFolderWavefront);
        end
    end

    fileFolderWavefront_net = [fileFolderIn '/Processed Cal Data_' expID '/wavefront_all_net/'];
    if SaveCM || SaveFunctionFigures || SaveScriptFigures
        if ~exist(fileFolderWavefront_net, 'dir')
            mkdir(fileFolderWavefront_net);
        end
    end

    %% Load and pre-process data

    aliststruct = dir([fileFolderIn '/A_*']);
    alist = str2double(arrayfun(@(x) string(x.name(3:end)), aliststruct));
    alist = sort(alist)';
    numA = length(alist);
    slopeM = zeros(numA,zernNum);
    fitRMSE=zeros(numA,zernNum);
    fitR2=zeros(numA,zernNum);
    phi_IF_all = zeros(Sy,Sx,numA);
    phi_IF_all_full = zeros(Sy,Sx,numA);
    IM_PD = zeros(length(zernikeIndicesForFitting),numA);
    IM_PD_full = zeros(length(zernikeIndicesIncludingPistonAndTilts),numA);
    fullc_A = zeros(numA,length(zernikeIndicesIncludingPistonAndTilts));

    %% Loop Over Actuators

    for A = alist
        % fprintf(['A = ', num2str(A), '\n']) %Uncomment to display A
        % number in real time
        pname = [fileFolderIn '/A_' num2str(A)];
        pokeliststruct = dir([pname '/A_Poke_*']);
        voltagelist = arrayfun(@(x) string(x.name), pokeliststruct);
        [uniquevoltagelist,~,voltagelistmatch] = unique(round(str2double(string(regexp(voltagelist,'(?<=Poke_)[-.\d+]+(?=.{4})','match'))),5));
        numV = length(uniquevoltagelist);

        if ~exist('imstackcrop','var')
            objshift = zeros(numV-1,4,numA); %Num. V - 1, 4 , num. A
            zernshift = zeros(numV-1,4,numA);
            imstackcrop = zeros(Sx,Sy,numz,numV,numA);
            objectEstimatestack = zeros(Sx,Sy,numV,numA);
            coeffsPD = zeros(numV,zernNum,numA,'single');
            phiRe_all = zeros(Sx,Sy,numV,numA);
        end

        coeffsPD_A = zeros(numV,zernNum);
        objectEstimatestackA = zeros(Sy,Sx,numV);
        objshiftA = zeros(numV-1,4); % Num. V - 1, 4 , num. A
        zernshiftA = zeros(numV-1,4);

        %% Loop over Voltages
        for V = 1:numV

            fileName=['A_' num2str(A) '_' num2str(uniquevoltagelist(V))];
            fileNameA=['A_' num2str(A)];
            fileFolderOutA = [fileFolderIn '/Processed Cal Data_' expID '/' fileNameA '/'];
            fileFolderOut = [fileFolderIn '/Processed Cal Data_' expID '/' fileNameA '/' fileName '/'];

            if SaveCM || SaveFunctionFigures || SaveScriptFigures
                if ~exist(fileFolderOut, 'dir')
                    mkdir(fileFolderOut);
                end
                if ~exist(fileFolderOutA, 'dir')
                    mkdir(fileFolderOutA);
                end
            end

            Vlistidx = find(voltagelistmatch==V);
            imstack=tiffreadVolume(strcat(pname,"/",voltagelist(Vlistidx)))-bgval;

            [cxy,~,~]=size(imstack);
            cx=cxy/2;
            cy=cx;

            imstackcrop(:,:,:,V,A) = imstack(cy-(HSy-1):cy+HSy,cx-(HSx-1):cx+HSx,:);

            if CMSH==true
                ImIn = fliplr(imstackcrop(:,:,:,V,A));
            else
                ImIn = imstackcrop(:,:,:,V,A);
            end

            phi_deltas = zeros(Sx, Sy, numz, 'single');
            RI = 1.33;
            phiunit = calc_defocusunit(Sx, pixelSize, lambda, NA, RI); %rad/um defocus?
            for i = 1:numz
                phi_deltas(:,:,i) = (ds(i)*phiunit);
            end

            %% Estimate Wavefront
            ii=ii+1;
            talgo=tic;
            [cEstimate, objectEstimate, ~] = reconstructZernikeAberrations(ImIn, zernikeIndicesForFitting, phi_deltas, ...
                gamma, 0, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU, 1, 0, fileFolderOut, SaveFunctionFigures);
            RunTime(ii)=toc(talgo);

            objectEstimatestack(:,:,V,A) = objectEstimate;
            objectEstimatestackA(:,:,V) = objectEstimate;

            %% Make Output Figures
            Vchar = char(voltagelist(Vlistidx));
            if SaveScriptFigures==true
                imwrite(rescale(objectEstimate),[fileFolderOut '/wavefront_a_' num2str(A) ' ' Vchar(8:end-4) ' ObjEst.png'])
                IFfig=figure();
                z3 = double(gather(zernikeObject.getAberration(zernikeIndicesForFitting, cEstimate) .* idx));
                phi_IF2 = 2*pi*z3/lambda; %um -> rad
                showWavefront = nan(Sx,Sy);
                showWavefront(idx) = phi_IF2(idx);
                pcolor(xi,xi,showWavefront);
                shading interp;
                axis square;
                xlim([minx maxx])
                ylim([miny maxy])
                clim([-5 5])
                axis off
                colormap parula
                Fc = colorbar;
                xlabel(Fc,'rad/V');
                xlim([minx-5 maxx+5])
                ylim([miny-5 maxy+5])
                title(['{\phi}(r,\theta): PhaseDiv wavefront_a_' num2str(A) ' ' Vchar(8:end-4) 'V']);
                drawnow
                exportgraphics(IFfig,[fileFolderOut '/wavefront_a_' num2str(A) ' ' Vchar(8:end-4) 'V.png'],'ContentType','vector')
                saveas(IFfig,[fileFolderOut '/wavefront_a_' num2str(A) ' ' Vchar(8:end-4) 'V.fig'])
                close(IFfig)
            end

            %% Calculate Tip/Tilt

            if V>1
                [tempobjectshift,~] = dftregistration(fft2(objectEstimatestackA(:,:,V-1)),fft2(objectEstimatestackA(:,:,V)),100);
                tempobjectshift(3:4) = tempobjectshift(3:4).*pixelSize;
                objshift(V,:,A) = tempobjectshift;
                objshiftA(V-1,:) = tempobjectshift;
                tempobjectshift(3:4) = tempobjectshift(3:4)./0.8352;
                zernshift(V,:,A) = tempobjectshift;
                zernshiftA(V-1,:) = tempobjectshift;
            end

            coeffsPD(V,:,A) = cEstimate;
            phi_c = 2*pi*cEstimate/lambda;
            phiRe = double(gather(zernikeObject.getAberration(zernikeIndicesForFitting, phi_c) .* idx));
            phiRe_all(:,:,V,A) = phiRe;
            coeffsPD_A(V,:) = cEstimate;
        end

        %% Compute Influence Function Slopes
        slopeM_A=zeros(1,zernNum);
        if SaveScriptFigures==true
            slopefig=figure();

            slopeax=tiledlayout('flow');
            slopeax.Padding='compact';
            slopeax.TileSpacing='loose';
        end

        for c=1:zernNum

            [xData, yData] = prepareCurveData(uniquevoltagelist,coeffsPD_A(:,c));
            [fitresult, gof] = fit( xData, yData, ft );
            fitMvals=coeffvalues(fitresult);
            fitRMSE(A,c)=gof.rmse;
            fitR2(A,c)=gof.rsquare;
            slopeM(A,c)=fitMvals(1);
            slopeM_A(c)=fitMvals(1);

            if SaveScriptFigures==true
                nexttile(c)
                hold on
                xlim([min(uniquevoltagelist) max(uniquevoltagelist)])
                title(['Zernike ' num2str(zernikeIndicesForFitting(c))])
                scatter(xData,yData,'filled','MarkerFaceColor','b')
                hline=refline(coeffvalues(fitresult));
                hline.LineWidth=1.5;
                hline.Color='c';
            end
        end

        if SaveScriptFigures==true
            title(slopeax,['Influence Function Slopes - Actuator ' num2str(A)],'FontSize',24)
            xlabel(slopeax,'Input Amplitude (Voltage)','FontSize',18)
            ylabel(slopeax,'Measured Coefficient Amplitude','FontSize',18)
            drawnow
            exportgraphics(slopefig,[fileFolderOutA '/Zernike Slopes ' fileNameA '.png'])
            saveas(slopefig,[fileFolderOutA '/Zernike Slopes ' fileNameA '.fig'])
            close(slopefig)
        end

        dx = uniquevoltagelist(end)-uniquevoltagelist(1);
        dc = slopeM_A;
        fullc=[mean(zernshiftA(:,3:4),1,'omitnan')./uniquetol(diff(uniquevoltagelist)) dc];
        fullc_A(A,:)=fullc;
        z = double(gather(zernikeObject.getAberration(zernikeIndicesForFitting, dc) .* idx));
        phi_IF = 2*pi*z/lambda;
        phi_IF_all(:,:,A) = phi_IF;

        z2 = double(gather(zernikeObject.getAberration(zernikeIndicesIncludingPistonAndTilts, fullc) .* idx));
        phi_IF2 = 2*pi*z2/lambda;
        phi_IF_all_full(:,:,A) = phi_IF2;

        if SaveScriptFigures==true
            IFfig=figure();
            showWavefront = nan(Sx,Sy);
            showWavefront(idx) = phi_IF(idx);
            pcolor(xi,xi,showWavefront);
            shading interp;
            axis square;
            axis off
            Fc = colorbar;
            xlabel(Fc,'rad/V');
            xlim([minx-5 maxx+5])
            ylim([miny-5 maxy+5])
            title(['{\phi}(r,\theta): PhaseDiv - ' num2str(length(uniquevoltagelist)) '-point fit']);
            colormap parula
            drawnow
            exportgraphics(IFfig,[fileFolderOutA '/wavefront_a_' num2str(A) ' ' num2str(numV) '-Point Slope.png'])
            saveas(IFfig,[fileFolderOutA '/wavefront_a_' num2str(A) ' ' num2str(numV) '-Point Slope.fig'])
            close(IFfig)

            IFfig=figure();
            showWavefront = nan(Sx,Sy);
            showWavefront(idx) = phi_IF2(idx);
            pcolor(xi,xi,showWavefront);
            shading interp;
            axis square;
            axis off
            Fc = colorbar;
            xlabel(Fc,'rad/V');
            xlim([minx-5 maxx+5])
            ylim([miny-5 maxy+5])
            title(['{\phi}(r,\theta): PhaseDiv with Tip/Tilt - ' num2str(length(uniquevoltagelist)) '-point fit']);
            colormap parula
            drawnow
            exportgraphics(IFfig,[fileFolderOutA '/wavefront_a_' num2str(A) ' ' num2str(numV) '-Point Slope with Tip-Tilt.png'])
            saveas(IFfig,[fileFolderOutA '/wavefront_a_' num2str(A) ' ' num2str(numV) '-Point Slope with Tip-Tilt.fig'])
            close(IFfig)

        end

        %% Compute and Save IM/CM

        if SaveCM==true
            if CMSH==true
                IM_PD(:,A) = coeffs_rot(dc, zernikeIndicesForFitting, -rotAng);
                IM_PD_full(:,A) = coeffs_rot(fullc, zernikeIndicesIncludingPistonAndTilts, -rotAng);
            else
                IM_PD(:,A) = dc;
                IM_PD_full(:,A) = fullc;
            end
        end
        % NOTE: taking this out of the loop would avoid recalculation of Zernikes!
    end



    if SaveCM==true

        if CMSH==true
            fileFolderOut = [fileFolderIn '/Processed Cal Data_' expID '/Cal Results/'];
        else
            fileFolderOut = [fileFolderIn '/Processed Cal Data_' expID '/Cal Results/nonrot'];
        end

        if ~exist(fileFolderOut, 'dir')
            mkdir(fileFolderOut)
        end


        CM_PD = pinv(IM_PD, invthresh);
        CM_PD_full = pinv(IM_PD_full, 0.005);

        if SaveScriptFigures==true
            PlotMirror(fileFolderOut,idx,phi_IF_all_full,xpup,ypup,xi,CM_PD_full,IM_PD_full')
        end

        dlmwrite([fileFolderOut 'IM_PD_full.dat'], IM_PD_full, 'delimiter', '\t', 'precision', '%.6f');
        dlmwrite([fileFolderOut 'CM_PD_full.dat'], CM_PD_full, 'delimiter', '\t', 'precision', '%.6f');

        if SaveMAT==true
        save([fileFolderOut '/CalResults.mat'],"CM_PD","CM_PD_full","IM_PD_full","IM_PD","phi_IF_all_full","fitRMSE","fitR2","slopeM","idx","coeffsPD","imstackcrop","phiRe_all","phi_IF_all");
        end
    end

    caltime(folderIndex) = toc(tStart);
    time_per_a(folderIndex) = caltime(folderIndex) ./ numA;
    time_per_v(folderIndex) = time_per_a(folderIndex) ./ numV;
end

fprintf(['\nRuntimes:\n', num2str(caltime), '\n'])
fprintf(['\nAverage time per A:\n', num2str(time_per_a), '\n'])
fprintf(['\nAverage time per V:\n', num2str(time_per_v), '\n'])
fprintf(['\nAverage time per V:\n', num2str(mean(caltime)), '\n'])
