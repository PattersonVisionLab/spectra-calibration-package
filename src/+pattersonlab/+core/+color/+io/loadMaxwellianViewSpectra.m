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

    if nargin < 5
        ndf = 0;
    end

    import pattersonlab.core.color.*;

    R = GammaRampMeasurement('660nm', folderName, values, beamDiameter, ...
        calibrationDate, 'MinWL', 620, 'MaxWL', 682, 'NDF', ndf, varargin{:});
    G = GammaRampMeasurement("530nm", folderName, values, beamDiameter, ...
        calibrationDate, 'MinWL', 525, 'MaxWL', 560, 'NDF', ndf, varargin{:});
    B = GammaRampMeasurement("420nm", folderName, values, beamDiameter, ...
        calibrationDate, 'MinWL', 410, 'MaxWL', 445, 'NDF', ndf, varargin{:});

    LED = LedCalibration([R, G, B]);
end