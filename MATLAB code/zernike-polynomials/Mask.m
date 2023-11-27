classdef Mask

    properties
        values double
        nGrid = 128
        mode = 'exact'  % 'exact' or 'FFT'
                        % 'exact':
                        % zero point is exactly in the middle,
                        % i.e. for even grid sizes between the two middle pixels
                        % 'FFT':
                        % zero point is the middle pixel for odd grid sizes
                        % and the pixel to the lower right of the exact centre for even grid sizes
                        % (this is also the pixel which MATLAB takes as zero)
        maxRadius = 1
        rotationAngle {mustBeInFullRadialRange(rotationAngle)} = 0 % in rad, clockwise rotation (mathematically negative!)
        radius
        theta
    end

    methods
        function obj = Mask(par)
            if nargin > 0
                obj = setInputParameters('Mask', obj, par);
            end

            % Calculate mask values
            N = [obj.nGrid, obj.nGrid];
            spacing = 2*obj.maxRadius / obj.nGrid;
            s = [spacing, spacing];
            if strcmp(obj.mode,'exact')==1
                x = linspace( -obj.maxRadius, obj.maxRadius, obj.nGrid );
                y = linspace( -obj.maxRadius, obj.maxRadius, obj.nGrid );
            elseif strcmp(obj.mode,'FFT')==1
                x = ( floor( -N(1)/2 + 0.5) : floor( N(1)/2 - 0.5) ) * s(1);
                y = ( floor( -N(2)/2 + 0.5) : floor( N(2)/2 - 0.5) ) * s(2);
            end

            [X, Y] = meshgrid(y,x);
            obj.radius = sqrt(X.^2+Y.^2);
            [~, angle] = getPolarCoordinates(obj);
            obj.theta = mod(angle + obj.rotationAngle + pi, 2*pi) - pi;
            
            obj.values = (obj.radius <= 1);
        end

        function [normalizedRadius, angle] = getPolarCoordinates(obj)
            % Normalized radius and theta
            x = linspace(-1,1,obj.nGrid);
            [X,Y] = meshgrid(x,x);
            [angle,normalizedRadius] = cart2pol(X,Y);
        end

        function plot(obj)
            imagesc(obj.values)
            set(gca,'YDir','normal')
            axis equal; axis tight;
        end
    end
end