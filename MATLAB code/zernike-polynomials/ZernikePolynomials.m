classdef ZernikePolynomials < handle

    properties %(Access=private)
        polynomials % ordering based on specified type (Noll or Ansi)
        mask Mask
        indexing {mustBeMember(indexing,{'NOLL','ANSI'})} = 'NOLL'
        normalization logical = 1
        gpu logical = 0
    end

    methods (Static)
        function obj = getInstance(mask, indexing, normalization, gpu)
            persistent instance
            if isempty(instance) || instance.useGpu~=gpu || ~isequal(size(instance.getMask), size(mask.values)) || ...
                    ~all(instance.getMask == mask.values,'all') || ~all(instance.getMaskTheta == mask.theta,'all')
                switch nargin
                    case 1
                        instance = ZernikePolynomials(mask);
                    case 2
                        instance = ZernikePolynomials(mask, indexing);
                    case 3
                        instance = ZernikePolynomials(mask, indexing, normalization);
                    case 4
                        instance = ZernikePolynomials(mask, indexing, normalization, gpu);
                end
            end
            obj = instance;
        end
    end

    methods (Access=private)
        function obj = ZernikePolynomials(mask, indexing, normalization, gpu)
            obj.mask = mask;
            [m,n] = size(mask.values);
            obj.polynomials = NaN(m,n,0);
            if nargin > 1
                obj.indexing = indexing;
            end
            if nargin > 2
                obj.normalization = normalization;
            end
            if nargin > 3
                obj.gpu = gpu;
            end
        end
    end

    methods (Access=public)
        function polynomials = getPolynomials(obj, linearIndices)
            if nargin < 2
                polynomials = obj.polynomials;
            else
                if obj.indexing == 'ANSI'
                    linearIndices = linearIndices + 1;
                    offset = 0;
                else
                    offset = 1;
                end
                maxIndexRequired = max(linearIndices);
                numberAlreadyCalculated = obj.getNumberCalculatedPolynomials;
                if (maxIndexRequired > numberAlreadyCalculated)
                    obj.polynomials = cat(3, obj.polynomials, ...
                        obj.calculatePolynomials((numberAlreadyCalculated+offset):maxIndexRequired));
                end
                polynomials = obj.polynomials(:,:,linearIndices);
            end
        end

        function mask = getMask(obj)
            mask = obj.mask.values;
        end

        function theta = getMaskTheta(obj)
            theta = obj.mask.theta;
        end

        function gpu = useGpu(obj)
            gpu = obj.gpu;
        end

        function numberCalculatedPolynomials = getNumberCalculatedPolynomials(obj)
            numberCalculatedPolynomials = size(obj.polynomials, 3);
        end

        function aberration = getAberration(obj, linearIndices, coefficients)
            assert(numel(linearIndices) == numel(coefficients))

            requiredPolynomials = obj.getPolynomials(linearIndices);
            % Weighted sum of Zernike polynomials
            weightedPolynomials = multiplyPagewise(requiredPolynomials, coefficients);
            aberration = sum( weightedPolynomials,3 );
        end

        function polynomials = calculatePolynomials(obj, linearIndices)
            [n,m] = obj.getZernikeIndices(linearIndices, obj.indexing);
            polynomials = ZernikePolynomials.calculateZernikePolynomials( n, m, obj.mask, obj.normalization, obj.gpu );
        end

        function plot(obj)
            for j = 1: obj.getNumberCalculatedPolynomials
                figure
                imagesc(obj.polynomials(:,:,j))
                colorbar; set(gca,'YDir','normal')
            end
        end
        
    end


    methods (Static)

        function [n,m] = getZernikeIndices(j, indexing)
            if nargin == 1
                indexing = 'NOLL';
            end
            switch indexing
                case 'NOLL'
                    n = floor(sqrt(2.*j) - 1/2);
                    s = mod(n,2);
                    m = (1-2.*mod(j,2)) .* (2 .* floor((1+2.*j-n.*(n+1)+s) ./ 4) - s);
                case 'ANSI'
                    d = sqrt(9+ 8*j);
                    n = ceil((d-3)/2);
                    m = 2*j - n.*(n+2);
            end
        end

        function polynomials = calculateZernikePolynomials(n, m, mask, doNormalize, gpu)
            % See https://en.wikipedia.org/wiki/Zernike_polynomials for
            % definition of Zernike polynomials and indexing conventions

            %% Check input
            if nargin<4
                doNormalize = 1;
            end

            assert( all(n>=abs(m)) && all(abs(m)>=0), 'Check input! It must hold that n>=|m|>=0.')

            %% Vectorize
            sizeMask = size(mask.values);
            rho = mask.radius(:);
            theta = mask.theta(:);
            mAbs = abs(m);

            %% Calculate polynomials

            [requiredPowersRho, powerIdx] = ZernikePolynomials.precalculateRequiredPowersOfRho(rho, n, m);

            polynomials = NaN(sizeMask(1),sizeMask(2),numel(n));
            for j = 1:numel(n)

                A = (n(j)+mAbs(j))/2;
                D = (n(j)-mAbs(j))/2;

                %% Radiual polynomials (sum from k = 0 to (n-m)/2)

                maxFactorial = max([n(j),A]);
                cumProds = cumprod([1,1:maxFactorial]);

                % (-1)^k (n-k)!
                nominator = (-1).^(0:D) .*  cumProds(n(j)-(0:D)+1);

                % k! * ((n+m)/2-k)! * ((n-m)/2-k)!
                F = cumProds((0:D)+1);
                denominator =  F .* cumProds(A-(0:D)+1) .* fliplr(F);

                % r^(n-2k)
                powers = n(j)-2*(0:D);
                powersRho = requiredPowersRho(:,powerIdx(powers+1));

                Rnm = sum(nominator./denominator .* powersRho,2);

                %% Zernike polynomials
                if m(j)==0
                    Z = Rnm;
                elseif m(j)>0 % 'even' Zernike polynomials
                    Z = Rnm .* cos(theta*mAbs(j));
                else % 'odd' Zernike polynomials
                    Z = Rnm .* sin(theta*mAbs(j));
                end

                %% Normalization
                if doNormalize
                    if m(j)==0
                        Z = Z *sqrt(n(j)+1);
                    else
                        Z = Z * sqrt(2*(n(j)+1));
                    end
                end

                polynomials(:,:,j) = reshape(Z,sizeMask).*mask.values;
            end

            if gpu
                polynomials = gpuArray(single(polynomials));
            end
        end

        function [requiredPowersRho, powerIdx] = precalculateRequiredPowersOfRho(rho, n, m)
            isEven = mod(n,2);
            if all(isEven) || all(~isEven)
                requiredPowers = min(abs(m)):2:max(n);
            else
                requiredPowers = min(abs(m)):1:max(n);
            end
            if requiredPowers(1)==0 % faster version if power 0 is required
                requiredPowersRho = arrayfun(@(p) rho.^p, requiredPowers(2:end), 'UniformOutput', false);
                requiredPowersRho = cat(2, requiredPowersRho{:});
                requiredPowersRho = [ones(length(rho),1) requiredPowersRho];
            else
                requiredPowersRho = arrayfun(@(p) rho.^p, requiredPowers, 'UniformOutput',false);
                requiredPowersRho = cat(2, requiredPowersRho{:});
            end
            powerIdx = NaN(max(n)+1,1);
            powerIdx(requiredPowers+1) = 1:length(requiredPowers);
        end
    end
end
