function [imgs, PSFs, waveFronts] = createPhaseDiversityImages(img0, p, coeffs, pixelSize, lambda, NA, nType, SNR, rotAng, flagGPU)
%
% Simulate phase diversity images based on ground truth image and phases
%
% Output
%   imgs: phase diversity images, third dimension: N, the number of images
%   waveFronts: wavefront images (in um unit) corresponding to input coeffs
%           The first image is the base or offset coefficients
% Input
%   img0: input ground truth image
%   p: a vector of single indexes(OSA/ANSI convention) for Zernike components,
%   coeffs : a matrix (unknown and known phases) of Zernike coefficients; the
%   first dimension corresponding to p; unit: um
%   the second dimension should be the image number N; first row is the
%   base or offset coefficients
%   pixelSize: pixel size of the intensity image, unit: um
%   lambda: wavelength, unit: um
%   NA: numerical aperture of objective
%   nType: noise type
%       1)'none': noise free
%       2)'gaussian': Gaussian noise
%       3)'poisson': Poisson noise
%       3)'both': Gaussion and Poisson noise
%   SNR: SNR
%   rotAng: rotation angle between wavefront (sensor) and diveristy
%       images, (unit: degree)

% By: Min Guo
% Jan 30, 2020

% Choice for SNR
SNRchoice = 2; % 1: define SNR based on the signal of ground truth
               % 2: define SNR based on the signal of aberrated image (1st image)
% if not define SNR (to be compatible to old version)
if nargin == 7
    SNR = 10;
    rotAng = 0;
elseif nargin == 8
    rotAng = 0;
end

[Sx, Sy] = size(img0);
imgNum = size(coeffs,1);
if Sx ~= Sy
    error('createImages: the x size of input images should be same as the y size');
end

if ~strcmp(nType, 'none')
    % desired signal level
    G = 10; % background white noise (Gaussian noise)
    S = addnoise_snr2signal(SNR, nType, G);

    % signal defined based on ground truth image
    if SNRchoice == 1
        sType = 3; % defined as the intensity average of the pixels > threshold.
        threPerc = 0.1; % threshold value
        S0 = addnoise_getSignal(img0,sType, threPerc);
    end
end

% Zernike coefficients: convert lengh unit(um) to phase unit(pi)
length2phase = 2*pi/lambda;
coeffs = length2phase * coeffs;

% Define the pupil coordinates (Polar coordinate system)
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
theta = theta + deg2rad(rotAng);
r0 = r(idx);
theta0 = theta(idx);

parMask.nGrid = Sx;
parMask.rotationAngle = deg2rad(rotAng);
parMask.maxRadius = calculatePupilRadius(Sx, pixelSize, lambda, NA);
pupilMask = Mask(parMask);
idx = logical(pupilMask.values);
zernikeObject = ZernikePolynomials.getInstance(pupilMask,'ANSI',0,flagGPU);

pupilMask = zeros(Sx, Sy,'single');
pupilMask(idx) = 1;
imgFFT = fft2(img0);
imgs = zeros(Sx, Sy, imgNum,'single');
PSFs = zeros(Sx, Sy, imgNum,'single');
waveFronts = zeros(Sx, Sy, imgNum,'single');
for i = 1:imgNum
    if i == 1
        c_phi = coeffs(1,:);
    else
        c_phi = coeffs(1,:) + coeffs(i,:);
    end
    phi = single(gather(zernikeObject.getAberration(p, c_phi) .* idx)); % in phase unit: pi
    pupilFun = pupilMask.*exp(1i*phi);
    prf = fftshift(ifft2(ifftshift(pupilFun)));
    PSF = abs(prf).^2;
    PSF = PSF/sum(PSF(:)); % sum normalized to 1: June 19, 2019
    OTF = fft2(ifftshift(PSF));
    imgBlurFFT = imgFFT.* OTF;
    imgBlur = real(ifft2(imgBlurFFT));
    if i == 1
        % signal defined based on aberrated image (1st image)
        if SNRchoice == 2
            sType = 3; % defined as the intensity average of the pixels > threshold
            threPerc = 0.1; % threshold value
            S0 = addnoise_getSignal(imgBlur,sType, threPerc);
        end
    end
    if strcmp(nType, 'none')
        img = imgBlur; % noise free
    else % scale all images by same value and add noise
        % image scaling ratio
        scaleR = S/S0;
        img = addnoiseScale(imgBlur, nType, scaleR, G);
    end

    imgs(:,:,i) = img;
    PSFs(:,:,i) = PSF;
    waveFronts(:,:,i) = phi;
end

% convert phase unit(pi) to lengh unit(um)
waveFronts = waveFronts/length2phase;
% remove offset
if imgNum > 1
    for i = 2: imgNum
        waveFronts(:,:,i) = waveFronts(:,:,i) - waveFronts(:,:,1);
    end
end

end
