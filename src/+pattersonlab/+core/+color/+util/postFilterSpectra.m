function spectra = postFilterSpectra(lightSpectra, filterSpectra)
% POSTFILTERSPECTRA
%
% Description:
%   Attenuates the spectral distribution of a light source by a filter
%
% Syntax:
%   spectra = calibration.util.postFilterSpectra(lightSpectra, filterSpectra)
%
% Inputs:
%   lightSpectra       Nx2 matrix
%       The spectral distribution of the light source (nm, normalized)
%   filterSpectra      Nx2 matrix or string
%       The fraction of light transmitted by a filter or the name of a
%       filter in the "data" directory.
%
% See also:
%   dichroicFilter, calibration.util.validateNDF
% --------------------------------------------------------------------------


    filterSpectra = convertCharsToStrings(filterSpectra);
    if isstring(filterSpectra)
        [~, filterSpectra] = pattersonlab.core.color.util.validateNDF(filterSpectra);
    end

    if any(filterSpectra(:, 2) > 1)
        % Convert from percent to fraction
        filterSpectra(:, 2) = filterSpectra(:, 2) / 100;
    end

    % Interpolate filter to match wavelengths of light source
    filt = interp1(filterSpectra(:,1), filterSpectra(:,2), lightSpectra(:,1));

    % Attenuate light source spectra by the fitler
    spectra = lightSpectra(:, 2) .* filt;

    % Remove NaNs
    spectra(isnan(spectra)) = 0;