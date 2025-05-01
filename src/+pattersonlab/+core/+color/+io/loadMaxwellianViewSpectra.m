function [LED, R, G, B] = loadMaxwellianViewSpectra(folderName, values, calibrationDate, beamDiameter, ndf, varargin)
% LOADMAXWELLIANVIEWSPECTRA
%
% Description:
%   Load LED calibration objects with default parameters for 1P AOSLO system
%
% Syntax:
%   [LED, R, G, B] = pattersonlab.core.color.loadMaxwellianViewSpectra(...
%       folderName, values, calibrationDate, beamDiameter, ndf, varargin)
%
% Additional key/value inputs are passed to GammaRampMeasurement.
%
% See also:
%   GammaRampMeasurement, LedCalibration
% --------------------------------------------------------------------------

    % arguments
    %     folderName      {mustBeFolder}
    %     values     (1,:) {mustBeVector}
    %     calibrationDate (1,1) string
    %     beamDiameter (1,1) double 
    %     ndf (1,1) double {mustBeNonnegative} = 0
    %     opts.RedMinWL = 620
    %     opts.RedMaxWL = 682
    %     opts.GreenMinWL = 525
    %     opts.GreenMaxWL = 560
    %     opts.BlueMinWL = 410
    %     opts.BlueMaxWL = 445
    % end


    import pattersonlab.core.color.*;

    R = GammaRampMeasurement('660nm', folderName, values, beamDiameter, ...
        calibrationDate, 'MinWL', 620, 'MaxWL', 682, 'NDF', ndf, varargin{:});
    G = GammaRampMeasurement("530nm", folderName, values, beamDiameter, ...
        calibrationDate, 'MinWL', 525, 'MaxWL', 560, 'NDF', ndf, varargin{:});
    B = GammaRampMeasurement("420nm", folderName, values, beamDiameter, ...
        calibrationDate, 'MinWL', 410, 'MaxWL', 445, 'NDF', ndf, varargin{:});

    LED = LedCalibration([R, G, B]);
end