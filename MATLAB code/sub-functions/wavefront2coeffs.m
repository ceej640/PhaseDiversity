function coeffs = wavefront2coeffs(waveFront, p, pupil, conType)
% wavefront2coeffs is to represent the wavefront by a series of Zernike coefficents with a given single-index p 
% Input
%   waveFront: an wavefront vector to be decomposited to zernike representations.
%   p: a vector of single indexes for Zernike components,
%       elements should be positive integers(>=1)
%   r: a vector of numbers between 0 and 1
%   theta: a vector of angles (rad), has same length with r
%   conType: single-index type, default: ANSI
% 
% Output
%   coeffs: a vector of zernike coefficients

if(nargin == 3)
    conType = 'ANSI';
end

zernike = ZernikePolynomials.getInstance(pupil,conType,0,0);
zernikeBasis = zernike.getPolynomials(p);
idx = logical(pupil.values);
zernikeBasis = reshape(zernikeBasis, [], size(zernikeBasis, 3));
zernikeBasis = zernikeBasis(idx(:),:);

coeffs = pinv(zernikeBasis)*waveFront;
coeffs = coeffs';
end
