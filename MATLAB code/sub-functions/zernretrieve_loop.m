function [cEstimate, imgEstimate, J, itTotal, cInter] = zernretrieve_loop(imgDs, ...
    r0, theta0, idx, p, waveFront_deltas, c0, zernPolynomials, Z, Rc, ...
    penalChioce, gamma, alpha, itLimit, stopChoice, tolValue, GPUflag, ~, pupilMaskObj)

% retrieve the zernike coefficients of aberrated phase based on phase
% diversity images;
% update formula(recurrence):
%                       c_i+1 = c_i + c_delta
%                       c_delta = -H^(-1)g
% Output
%   cEstimate: retrieved Zernike coefficients(um)
%   imgEstimate: object estimate
%   J: cost function value at each iteration
%   itTotal: number of iterations required
%   cInter: Zernike coefficients at each iteration
% 
% Input
%   imgDs: FT of phase diversity images, third dimension: N, the number of images
%   r0,theta0,idx: define the pupil aperture of the wavefront
%   1)r0: a 2D matrix of numbers between 0 and 1
%   2)theta0: a 2D matrix of angles (rad), has same size with r0
%   3)idx:
%       elements should be positive integers(>=1)
%   p: a vector of single-index Zernike coefficients,
%   waveFront_deltas: wavefronts corresponding to phase diversity images
%   (i.e. known/applied aberrations)
%   Rc: Rc matrix (unused when alpha = 0, as typical)
%   penalChioce: penalty term to image,
%       1 for L2 norm of image and 2 for L2 norm of gradient. Typical = 1
%   gamma: parameter for the regularization to image
%   alpha: not needed, set to 0.
%   itLimit: maximum allowed iteration number
%   stopChoice: iteration stopping criterion
%       0: no stopping criterion;
%       1: based on RMS of wavefront;
%       2: based on loss function;
%   tolValue = Typically = 0.001; % iteration stopping criterion for change
%   in objective function J.
%   GPUflag: GPU options, 0: CPU; 1: GPU;


% By: Min Guo
% Jan 29, 2020
% Modified by: Magdalena Schneider

r = r0(idx); % convert 2D matrix to vector with selected elements
theta = theta0(idx); % convert 2D matrix to vector with selected elements


[Sx, Sy, imgNum] = size(imgDs); % image size and numbers
nZernike = length(p);
idxNum = length(r);
J = zeros(itLimit+1, 1, 'single'); % cost function values (initial + all iterations)
cInter = zeros(itLimit, nZernike, 'single');% Zernike coefficients (all iterations)
Jrms = zeros(itLimit, 1, 'single');
% GPU Data transfering
if(GPUflag==1)

    pupilMask = zeros(Sx, Sy, 'single', 'gpuArray');

    cEstimate = gpuArray(c0);
    cDelta = zeros(size(c0),'single', 'gpuArray');

    g = zeros(nZernike,1, 'single', 'gpuArray');
    J = gpuArray(J);
    Jrms = gpuArray(Jrms);
    cInter = gpuArray(cInter);
else

    pupilMask = zeros(Sx, Sy, 'single');


    cEstimate = c0;
    cDelta = zeros(size(c0),'single');
    g = zeros(nZernike,1, 'single');
end
pupilMask(idx) = 1;

% calculate penalty term
switch penalChioce
    case 1
        penalTerm = gamma;
    case 2
        penalTerm = gamma* r0.*r0; % gamma*u^2
    otherwise
        error('zernretrieve_loop: wrong penalty choice');
end

% for calculating cost function
imgDsSq = sum(abs(imgDs).^2,3);

pairs = nchoosek(1:imgNum,2);
pairs = sortrows(fliplr(pairs));

zernikeObject = ZernikePolynomials.getInstance(pupilMaskObj,'ANSI',0,GPUflag);
Z = reshape(Z,Sx,Sx,1,nZernike);

for i = 1:itLimit+1
   
    % ================== %
    % Gradient: g
    % ================== %

    c = cEstimate;
    waveFront = zernikeObject.getAberration(p, c) .* idx;

    phis = waveFront + waveFront_deltas; % unknown (estimated) + known (phase diversity) aberrations
    Hks = pupilMask.*exp(1i*phis); % pupil function
    hks = ifft2(Hks); % prh
    sks = abs(hks).^2; % PSF
    Sks = fft2(sks); % FT of PSF

    % calculate Q
    Q = sum(abs(Sks).^2,3) + penalTerm;

    % % calculate cost function value: last iteration

    phaseTerm = alpha*c'*Rc*c;
    numeTerm = sum(conj(imgDs).*Sks,3);
    tTerm = imgDsSq - (abs(numeTerm).^2)./Q;
    J(i) = sum(tTerm(:)) + phaseTerm;

    % % % % check if terminate iteration
    %       0: no stopping criterion;
    %       1: based on RMS of wavefront;
    %       2: based on loss function;

    switch(stopChoice)
        case 1 % based on RMS of wavefront;
            waveFrontTemp = zernPolynomials*c;
            Jrms(i) = rms(waveFrontTemp);
            waveFrontTemp = zernPolynomials*cDelta;
            rmsDelta = rms(waveFrontTemp);
            if(i>=3)
                if(rmsDelta >= Jrms(i)||rmsDelta >= Jrms(i-1))
                    itTotal = i-2; % diverges, use second last estimate
                    break;
                end
                if(rmsDelta/Jrms(i)<tolValue)
                    itTotal = i-1; % converges, use last estimate
                    break;
                end
            end
        case 2 % based on RMS of wavefront;
            Jlast = J(i);
            if(i>=3)
                %               % compare to last two iterations to avoid fluctuations
                if((J(i-2)<Jlast)&&(J(i-1)<Jlast)) % diverges, use second last estimate
                    itTotal = i-2;
                    break;
                end
            end
            if(i>=2)
                Jdelta = abs(J(i-1) - Jlast)/abs(J(1)-Jlast);
                if(Jdelta<tolValue)
                    itTotal = i-1; % converges, use last estimate
                    break;
                end
            end
        otherwise % no stopping criterion;
    end
    % if reach maximum iteration
    if(i==itLimit+1)
        itTotal = itLimit;
        break;
    end


    % calculate F
    numeTerm = sum(conj(Sks).*imgDs,3);
    F = numeTerm./Q;

    % calculate Vk and g[phi]
    g_phi = zeros(idxNum,1, 'single');
    if(GPUflag==1)
        g_phi = gpuArray(g_phi);
    end
    Vk = conj(F).*imgDs - abs(F).^2.*Sks;
    Vk_iFT = ifft2(Vk);
    temp1 = hks.*real(Vk_iFT);
    temp2 = fft2(temp1);
    temp3 = imag(conj(Hks).*temp2);
    temp4 = sum(temp3,3);
    g_phi = g_phi - 2*temp4(idx);

    % calculate gradient: g
    g = (g_phi' * zernPolynomials)' + 0.5 * alpha * Rc * c;

    % ======================== %
    % Hessian matrix: Hmatrix
    % ======================== %
    % simplify to K=2: k=2,j=1

    DQks = imgDs ./ sqrt(Q);

    hks_conj = conj(hks);

    hks_j = hks(:,:,pairs(:,2)); % j-loop
    hks_k = hks(:,:,pairs(:,1)); % k-loop
    Hks_j = Hks(:,:,pairs(:,2)); % j-loop
    Hks_k = Hks(:,:,pairs(:,1)); % k-loop
    DQks_j = DQks(:,:,pairs(:,2)); % j-loop
    DQks_k = DQks(:,:,pairs(:,1)); % k-loop
    
    Hks_phi_n = Hks .* Z;
    Hks_ifft = ifft2(Hks_phi_n);
    TEMP = fft2( imag(hks_conj .* Hks_ifft) );
    TEMP_k = TEMP(:,:,pairs(:,1),:); % k-loop
    DQks_mult_TEMP_jk = DQks_j.*TEMP_k;
    TEMP_j = TEMP(:,:,pairs(:,2),:); % j-loop
    DQks_mult_TEMP_kj = DQks_k.*TEMP_j;

    U = DQks_mult_TEMP_jk - DQks_mult_TEMP_kj;

    temp3 = conj(Hks_j) .* fft2( hks_j .* ifft2(conj(DQks_k).*U) );
    temp6 = conj(Hks_k) .* fft2( hks_k .* ifft2(conj(DQks_j).*U) );
    HGN_phi_n = 4* imag(temp3-temp6);
    HGN_phi_ns = squeeze(sum(HGN_phi_n,3));

    % calculate Hessian Matrix elements
    HGN_phi_ns_reshaped = reshape(HGN_phi_ns,[],nZernike);
    HGN_phi_ns_relevant = HGN_phi_ns_reshaped(idx,:);
    Hmatrix = HGN_phi_ns_relevant' * zernPolynomials;

    Hmatrix = Hmatrix + alpha*Rc;
    cDelta = - Hmatrix\g;
    cEstimate = c + cDelta;
    cInter(i, :) = cEstimate;
end
imgEstimate = real(ifft2(F));
imgEstimate(imgEstimate<0) = 0;
cEstimate = cInter(itTotal, :);
if (GPUflag == 1)
    J = gather(J);
    cEstimate = gather(cEstimate);
    imgEstimate = gather(imgEstimate);
    cInter = gather(cInter);
end

