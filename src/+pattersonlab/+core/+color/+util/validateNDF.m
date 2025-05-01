function [tf, spectra] = validateNDF(ndfName)
% VALIDATENDF
%
% Description:
%   Validates the name of a filter file and loads the spectra if valid.
%
% Syntax:
%   [tf, spectra] = validateNDF(ndfName)
% --------------------------------------------------------------------------

    filterDir = fullfile(pattersonlab.core.color.util.getMainFolder(), 'data');

    fNames = getFolderFiles(filterDir);
    fNames = fNames(startsWith(fNames, "NE"));

    ndfName = erase(ndfName, ".txt");
    fNames = erase(fNames, ".txt");

    if ~ismember(ndfName, fNames)
        warning("NDF file named %s not found\n", ndfName);
        disp(fNames)
        tf = false;
    else
        tf = true;
    end

    if nargout == 2
        if tf
            spectra = dlmread(fullfile(filterDir, ndfName + ".txt"));
            % Convert from percent to fraction (assuming NDF)
            spectra(:, 2) = spectra(:, 2) / 100;
        else
            spectra = [];
        end
    end
