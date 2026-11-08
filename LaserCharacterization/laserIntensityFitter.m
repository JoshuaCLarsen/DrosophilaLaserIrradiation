% Author: Joshua Larsen
%
% Description: Generates a simulation of laser power over the region of interest,
%              approximating the axes in mm instead of pixels. The default
%              imageDistance and focalLength values match the manuscript setup.

function laserFit = laserIntensityFitter(imageDistance, focalLength)

arguments
% Default arguments were used in the manuscript; both lengths are in meters.
    imageDistance (1, 1) double {mustBeFinite, mustBePositive} = 1.0
    focalLength (1, 1) double {mustBeFinite, mustBePositive} = 35e-3
end

laserImagePath = "rawLaserIntensityPattern.png";
laserImageValues = imread(laserImagePath);
laserImageBlues = double(laserImageValues(:,:,3)) / 255;

%% Do some scaling to generate pixel size and array dims
pixelPitch = 1.55e-3; % Measured in mm
pixelPitchM = 1.55e-6; % Measured in m
[ny, nx] = size(laserImageBlues);
% x = (0:nx-1) * pixelPitch;
% y = (0:ny-1) * pixelPitch;

% We also need the ground sampling distance
gsd = imageDistance * pixelPitchM / focalLength * 1000;
x = (0:nx-1) * gsd;
y = (0:ny-1) * gsd;

%% Generate the Image and save it off
imagesc(x, y, 255 * laserImageBlues); colormap('jet');
colorbar;
axis image
xlabel('Horizontal Position (mm)');
ylabel('Veritcal Position (mm)');
c = colorbar('Ticks', [0:50:250]);
c.Label.String = 'Intensity (Digital Numbers)';
title('Laser Power Spatial Distribution');

% Save off the image for use in our paper
exportgraphics(gcf, "laserImage.pdf");

%% Generating the fitted version based on fitting performed for the paper
z = laserImageBlues;
[X, Y] = meshgrid(1:nx, 1:ny);
%fittedLaserProp = fit([X(:), Y(:)], 255 * laserImageBlues, 'poly2');
% Create the custom supergaussian function here
ft = fittype( ...
'A * exp(-((x - xOrigin) ./ wx) .^ 2 -((y - yOrigin) ./ wy) .^ 6)', ...
'independent',{'x','y'}, ...
'dependent', 'z',...
'coefficients',{'A', 'xOrigin','yOrigin','wx','wy'});

% Providing an initial guess to help with fitting.
AOrigin = 0.95;
xOrigin = nx/2;
yOrigin = ny/2;
wxOrigin = nx/5;
wyOrigin = ny/5;

start = [AOrigin xOrigin yOrigin wxOrigin wyOrigin];
% Actually run the fit now that everything has been set up.
[fitresult, gof] = fit([X(:), Y(:)], z(:), ft, 'StartPoint', start);

% We need to save the resulting values to a table that can be written out
% to a csv or xlsx file.
laserFit.A = 255 * fitresult.A;
laserFit.WX = fitresult.wx .* gsd;
laserFit.WY = fitresult.wy .* gsd;
laserFit.Y0 = fitresult.yOrigin .* gsd;
laserFit.X0 = fitresult.xOrigin .* gsd;
laserFit.NX = 2;
laserFit.NY = 6;
fittedLaser = struct2table(laserFit);
writetable(fittedLaser,"LaserFittingResults.csv");

% Generate the accompanying fit on pixels.
imagesc(x, y, 255 * abs(fitresult(X, Y)), [0 255]);
c = colorbar('Ticks', [0:50:250]);
c.Label.String = 'Intensity (Digital Numbers)';
axis image
xlabel('Horizontal Position (mm)');
ylabel('Veritcal Position (mm)');
title('Simulated Laser Intensity');
exportgraphics(gcf, "LaserFit.pdf");

% Generate the accompanying fit on pixels, then get the residuals.
imagesc(x, y, 255 * abs(fitresult(X, Y) - laserImageBlues), [0 255]);
c = colorbar('Ticks', [0:50:250]);
c.Label.String = 'Intensity (Digital Numbers)';
axis image
xlabel('Horizontal Position (mm)');
ylabel('Veritcal Position (mm)');
title('Difference between Laser Image and Simulated Intensity Pattern');
exportgraphics(gcf, "LaserResiduals.pdf");
end


function superGausDist = superGaus(x, xw, xOrigin, xn, y, yw, yOrigin, yn)
[X, Y] = meshgrid(x, y);
superGausDist = exp(-((X - xOrigin) / xw) .^ xn) .* exp(-((Y - yOrigin) / yw) .^ yn);
end

