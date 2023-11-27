function coeffs = coeffs_rot(coeffs0, p, rotAng)
% coeffs_rot is to convert Zernike coefficients by rotating wavefront
% Output
%   coeffs: a vector of Zernike coefficients
% Input
%   p: a vector of single indexes for Zernike components,
%   coeffs0: a vector of input Zernike coefficients corresponding to p
%   r: a vector of numbers between 0 and 1
%   rotAng: rotation angle, unit: degree

parMask.nGrid = 256;
parMask.rotationAngle = deg2rad(mod(rotAng,360));
parMask.maxRadius = 1;
pupilMask = Mask(parMask);
idx = logical(pupilMask.values);
zernikeObject = ZernikePolynomials.getInstance(pupilMask,'ANSI',0,0);
waveFrontVector = single(gather(zernikeObject.getAberration(p, coeffs0) .* idx));
waveFrontVector = flipud(waveFrontVector);

waveFrontVector = waveFrontVector(idx);

parMask.rotationAngle = deg2rad(mod(0,360));
pupilMaskUnrotated = Mask(parMask);
coeffs = wavefront2coeffs(waveFrontVector(:), p, pupilMaskUnrotated);
