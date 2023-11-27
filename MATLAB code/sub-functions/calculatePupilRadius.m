function pupilRadius = calculatePupilRadius(Sx, pixelsize, wavelength, NA)
%Calculate Pupil Radius
%
%output: Pupil radius
%
%input:
%
%Sx: Image size
%
%pixelSize: Image pixel size in um
%
%wavelength: Central emission wavelength in um
%
%NA: Numerical aperture
%
freMax = NA/wavelength; % length^-1
freSampling = 1/pixelsize; % length^-1
freSamplingPhase = Sx/freSampling;
pixelSizePhase = 1/freSamplingPhase; % 1/(Sx*pixelSize)
pupilRadius = round(Sx/2)*pixelSizePhase/freMax;
end