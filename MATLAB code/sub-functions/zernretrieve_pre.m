function [imgDs, zernikeBasis, Rc, waveFront_deltas] = zernretrieve_pre(imgs, r0, theta0, idx, p, c_delta, nFlag, gpuFlag, pupilMask)
% perform pre-calculations for wavefront estimation
%
% Output

%   imgDs: FT of phase diversity images, third dimension: N, the number of images
% 
%   zernikeBasis: basis of Zernike components
% 
%   Rc: unused when alpha = 0 as typical.
% 
%   waveFront_deltas: wavefronts corresponding to phase diversity images
%  
% Input
% 
% imgs: phase diversity images, third dimension: N, the number of images
% 
% r0,theta0,idx: define the pupil aperture of the wavefront
%   1)r0: pupil radius: a 2D matrix of numbers between 0 and 1
%   2)theta0: a 2D matrix of angles (rad), has same size with r0
%   3)idx:
%       elements should be positive integers(>=1) corresponding to indices
%       in the range of the pupil. 
% 
% p: a vector of single indexes for Zernike components,
% 
% c_delta: a vector or matrix (known aberration) of Zernike coefficients
%   corresponding to p (should be normalized to phase unit: pi); if matrix,
%   the second dimension should be N-1

% By: Min Guo
% Jan 29, 2020
% Modifications: July 1, 2020
%   modify input pupil parameters(r, theta) from 1D vectors to 2D matrix 

r = r0(idx); % convert 2D matrix to vector with selected elements
[Sx, Sy, imgNum] = size(imgs); % image numbers:
deltaNum = size(c_delta,1);
if (imgNum - deltaNum)~= 1
    error('zernretrieve_pre:NMlength','deltaNum should be: imgNum -1.')
end
idxNum = length(r);
opNum = length(p);

imgDs = single(fft2(imgs)); % FT of imgs

% Calculate Rc
Zernike = ZernikePolynomials.getInstance(pupilMask,'ANSI',0,gpuFlag);
zernikeBasis = single(Zernike.getPolynomials(p)); % Zernike polynomial basis

zernP2Dtemp1 = zeros(Sx, Sy, opNum, 'single');
zernP2Dtemp2 = zeros(Sx, Sy, opNum, 'single');
zernP2Dtemp1(:,2:Sy,:) = zernikeBasis(:,1:Sy-1,:);
zernP2Dshifted = zernP2Dtemp1 - zernikeBasis;
zernP2Dshifted = reshape(zernP2Dshifted, [], size(zernP2Dshifted, 3));
zernPhs = zernP2Dshifted(idx(:),:);
zernP2Dtemp2(2:Sx,:,:) = zernikeBasis(1:Sx-1,:,:);
zernP2Dshifted = zernP2Dtemp2 - zernikeBasis;
zernP2Dshifted = reshape(zernP2Dshifted, [], size(zernP2Dshifted, 3));
zernPvs = zernP2Dshifted(idx(:),:);

Rc = zernPhs' * zernPhs + zernPvs' * zernPvs;


%% Calculate wavefront aberrations

waveFront_deltas = zeros(Sx,Sy,imgNum, 'single'); % known aberration, e.g., defocus
for k = 2:imgNum
    waveFront_deltas(:,:,k) = single(gather(Zernike.getAberration(p, c_delta(k-1,:)) .* idx));
end

