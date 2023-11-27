function imageQualityMetrics = calculateImageQualityMetrics(image, binSize)
%
% Calculate image quality metrics
% Metric type: norm DCTS (discrete cosine transform)  - spectral domain 3: 0~r0
%
% Output:
%   imageQualityMetrics
% Input:
%   image: input image, same image size in x and y
%   binSize

if (nargin == 2) && (binSize ~= 1)
    flagBin = 1;
else
    flagBin = 0;
end

wavelength = 0.532;
pixelSize = 0.104;
NA = 1.2;

if flagBin
    image = imbin(image, binSize);
    pixelSize = pixelSize * binSize;
end

[Sx,~] = size(image);
maxFrequency = 2*NA / wavelength;
pixelSizeFrequencySpace = 1 / (pixelSize*Sx);
r0 = maxFrequency / pixelSizeFrequencySpace;
imgDCT = dct2(image);
imgNorm = norm(imgDCT);
imgABS = abs(imgDCT / imgNorm);

[X, Y] = meshgrid(1:Sx);
isWithinRadius = ((X+Y) <= r0);
isNonZeroImgDCT = (imgDCT ~= 0);
isIndexForSum = (isWithinRadius & isNonZeroImgDCT);
imgAbs = imgABS(isIndexForSum);
sumValue = sum(imgAbs .* log2(imgAbs), 'all');

imageQualityMetrics = - 2/r0^2 * sumValue;

end