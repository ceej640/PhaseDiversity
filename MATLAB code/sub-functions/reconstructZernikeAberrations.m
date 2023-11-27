function [cEstimate, imgEstimate, rePar] = reconstructZernikeAberrations(images, p, aberration_delta, gamma,...
    alpha, itLimit, zernSteps, pixelSize, lambda, NA, GPUflag, penalChoice, rotAng, savedir, plotFigures)
% reconstruct the zernike coefficients of aberrated phase
% based on phase diversity images and upboosting strategy
%
% Output
%   cEstimate: reconstructed zernike coefficients(um)
%   imgEstimate: estimated image
%   rePar: record of intermediate results
% Input
%   images: phase diversity images, third dimension: N, the number of images
%   p: Diversity aberration information. Input as a vector of single-index Zernike coefficients to be estimated,
%   aberration_delta: either a vector of Zernike coefficients or matrix of wavefront aberrations
%   corresponding to p (unit: um); if matrix, the second dimension should be N-1
%   gamma: parameter for the regularization to image
%   alpha: parameter for the regularization to phase
%   itLimit: maximum iteration number
%   zernSteps: coefficients in p to be updated step by step. Typically set
%   to a single value equal to the number of coefficients being estimated
%   pixelSize: pixel size of the intensity image, unit: um
%   lambda: wavelength, unit: um
%   NA: numerical aperture of objective
%   GPUflag: GPU options, 0: CPU; 1: GPU;
%   penalChioce: penalty term to image,
%       typically = 1 = penalize square norm
%   rotAng: rotation angle between wavefront (sensor) and diveristy
%       images, (unit: degree). Typically set to 0.

%%

[Sx, Sy, imgNum] = size(images);
if Sx ~= Sy
    error('recon_zern: the x size of input images should be same with the y size');
end

if length(p) < max(zernSteps(:))
    error('recon_zern: upSteps exceed the p index');
end

if size(aberration_delta,1) == numel(p) || size(aberration_delta,2) == numel(p) % vector of Zernike coefficients
    c_delta = aberration_delta;
    aberrationInputType = 'Zernike coefficients';
    deltaNum = size(c_delta,1);
    if (imgNum - deltaNum) ~= 1
        error('recon_zern:NMlength','deltaNum should be: imgNum -1.')
    end
else % matrices of wavefront aberrations
    waveFront_deltas = aberration_delta;
    aberrationInputType = 'Wavefront aberrations';
    deltaNum = size(waveFront_deltas,3);
    if imgNum~= deltaNum
        error('recon_zern:NMlength','deltaNum should be: imgNum.')
    end
end


%%

if nargin == 11
    penalChoice = 1; % default for Tikhonov norm
    rotAng = 0;
elseif nargin == 12
    rotAng = 0;
end

if nargin < 15
    plotFigures = 0;
end

% % % customize a few settings:
nFlag = 'none'; % Zernike norm option
% iteration stopping criterion
%       0: no stopping criterion;
%       1: based on RMS of wavefront; typically tolValue = 0.01;
%       2: based on loss function; typically tolValue = 0.001;
stopChoice = 2;
tolValue = 0.001; % tolerance value

% Define the pupil coordinates (Polar coordinate system)
parMask.nGrid = Sx;
parMask.rotationAngle = deg2rad(rotAng);
parMask.maxRadius = calculatePupilRadius(Sx, pixelSize, lambda, NA);
pupilMask = Mask(parMask);
theta0 = pupilMask.theta;
r0 = pupilMask.radius;
idx = logical(pupilMask.values);
zernikeObject = ZernikePolynomials.getInstance(pupilMask,'ANSI',0,GPUflag);

% Zernike coefficients: convert lengh unit(um) to phase unit(pi)
length2phase = 2*pi/lambda;
switch aberrationInputType
    case 'Zernike coefficients'
        c_delta = length2phase * c_delta;
        % Input diversity phases as Zernike coefficients
        [imgDs, zernikeBasis, Rc, waveFront_deltas] = zernretrieve_pre(images, ...
            r0, theta0, idx, p, c_delta, nFlag, GPUflag, pupilMask);
    case 'Wavefront aberrations'
        c_delta = zeros(imgNum-1, length(p), 'single');
        % Input diversity phases as wavefront images
        [imgDs, zernikeBasis, Rc, ~] = zernretrieve_pre(images, ...
            r0, theta0, idx, p, c_delta, nFlag, GPUflag, pupilMask);
end

% reshape
opNum = length(p);
zernPolynomials = NaN(sum(idx(:)), opNum);
for n = 1:opNum
    tmp = zernikeBasis(:,:,n);
    tmp = tmp(idx);
    zernPolynomials(:,n) = tmp;
end

% reconstruction
stepNum = length(zernSteps);
rePar.cEstimate = zeros(length(p),stepNum, 'single');
rePar.J = zeros(itLimit+1,stepNum, 'single');
rePar.itTotal = zeros(1,stepNum, 'single');
rePar.cInter = cell(1,stepNum);
cEstimate = zeros(zernSteps(1),1, 'single');
if GPUflag == 1
    imgDs = gpuArray(imgDs);
    zernPolynomials = gpuArray(zernPolynomials);
    Rc = gpuArray(Rc);
    waveFront_deltas = gpuArray(waveFront_deltas);
end


if plotFigures
    xi =1:Sx;

    rms_diversity = [];
    for k = 1:imgNum
        rms_diversity(k) = rms(rms(waveFront_deltas(:,:,k)));
    end

    [xx,yy] = ind2sub(size(idx),find(idx));

    hF = figure;
    ha = axes;
    hold on

    hFc = figure;
    hac = axes;
    hold on

    c_cat = [];

    hWF = figure;
    haxWF = tiledlayout('flow');
    haxWF.TileSpacing = 'tight';
    haxWF.Padding = 'tight';
    title(haxWF,'{\phi}(r,\theta): PhaseDiv - wavefronts by iteration');
end


% stepNum = 1 when graduated optimization is not used (as typical)
for iStep = 1:stepNum
    zernN = zernSteps(iStep);
    p0 = p(1:zernN);
    c0 = zeros(zernN,1, 'single');
    zernNLastIt = length(cEstimate);
    c0(1:zernNLastIt) = cEstimate;
    zernPolynomials0 = zernPolynomials(:,1:zernN);
    Z0 = zernikeBasis(:,:,1:zernN);
    Rc0 = Rc(1:zernN,1:zernN);
    [cEstimate, imgEstimate, J, itTotal, cInter] = zernretrieve_loop(imgDs,...
        r0, theta0, idx, p0, waveFront_deltas, c0, zernPolynomials0, Z0, Rc0,...
        penalChoice, gamma, alpha, itLimit, stopChoice, tolValue, GPUflag, nFlag, pupilMask);
    rePar.cEstimate(1:zernN,iStep) = cEstimate;
    rePar.J(:,iStep) = J;
    rePar.itTotal(iStep) = itTotal;
    rePar.cInter{iStep} = cInter;

    if plotFigures
        cIntercat = NaN(itTotal,max(zernSteps));

        cInter(~all(cInter,2),:) = [];
        cIntercat(1:size(cInter,1),1:zernSteps(iStep)) = cInter;

        xaxmin = sum(vertcat(rePar.itTotal)) - rePar.itTotal(iStep) + 1;
        xaxmax = xaxmin+length(J(J~=0))-1;
        plot(ha,xaxmin:xaxmax,J(J~=0))

        c_cat = vertcat(c_cat,cIntercat);

        rms_measured = [];
        for q = 1:size(cInter,1)
            showWavefront = nan(Sx,Sy);
            z = double(gather(zernikeObject.getAberration(p0, cInter(q,:)) .* idx));
            if any(z(:))
                rms_measured(q) = rms(rms(z(idx)));
                nexttile(haxWF)

                showWavefront(idx) = z(idx);

                pcolor(xi,xi,showWavefront), shading interp
                axis square;
                axis off;

                xlim([min(xx(:))-5 max(xx(:))+5])
                ylim([min(yy(:))-5 max(yy(:))+5])

                title(['it. ' num2str(itLimit*(iStep-1)+q)]);

                minz(q) = min(z(idx));
                maxz(q) = max(z(idx));
            end
        end
        minminz(iStep) = min(minz);
        maxmaxz(iStep) = max(maxz);
    end
end

% convert phase unit(pi) to lengh unit(um)
cEstimate = cEstimate/length2phase;
rePar.cEstimate = rePar.cEstimate/length2phase;


%% Plot figures

if plotFigures

    figure(hF)
    axes(ha)
    xlabel('Iteration Number')
    ylabel('J (Minimization Term')
    title('Minimization of J')
    shg
    drawnow
    exportgraphics(hFc,[savedir '/Minimization of J - Max. ' num2str(itLimit) 'Iterations.png'])
    saveas(hF,[savedir '/Minimization of J - Max. ' num2str(itLimit) 'Iterations.fig'])
    close(hF)


    figure(hFc)
    axes(hac)
    itaxis = 1:size(c_cat,1);
    plot(hac,itaxis,c_cat)
    xlabel('Iteration Number')
    ylabel('Zernike Coeff. Amplitude (Minimization Term')
    title('Convergence of Zernikes')

    drawnow
    exportgraphics(hFc,[savedir '/Convergence of Zernikes - Max. ' num2str(itLimit) 'Iterations.png'])
    saveas(hFc,[savedir '/Convergence of Zernikes - Max. ' num2str(itLimit) 'Iterations.fig'])
    close(hFc)

    figure(hWF)
    wMin_c = min(minminz);
    wMax_c = max(maxmaxz);
    for q = 1:(size(haxWF.Children,1)/2)
        nexttile(q)
        if wMax_c > wMin_c
            clim([wMin_c wMax_c])
        end
    end

    Fc = colorbar;
    Fc.Layout.Tile = 'east';
    Fc.Label.String = 'rad';
    Fc.Limits = [wMin_c wMax_c];
    drawnow
    exportgraphics(hWF,[savedir '/wavefronts by iteration - Max. ' num2str(itLimit) 'Iterations.png'])
    saveas(hWF,[savedir '/wavefronts by iteration - Max. ' num2str(itLimit) 'Iterations.fig'])
    close(hWF)
end

