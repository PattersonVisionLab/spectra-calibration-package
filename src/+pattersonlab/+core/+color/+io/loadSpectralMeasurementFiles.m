function [data, fileNames] = loadSpectralMeasurementFiles(fPath, token)
% LOADSPECTRALMEASUREMENTFILES
% 
% Syntax:
%   [data, fileNames] = loadSpectralMeasurementFiles(fPath, token)
% -------------------------------------------------------------------------

    arguments
        fPath           string      {mustBeFolder}
        token           string      = ""
    end

    fileNames = deblank(string(ls(fPath)));
    fileNames = fileNames(endsWith(fileNames, ".txt"));
    if token ~= ""
        fileNames = fileNames(contains(fileNames, token));
    end

    if all(isempty(fileNames)) || all(fileNames == "")
        error('No files found with token %s in %s', token, fPath);
    end

    % Sort them
    shortNames = cellfun(@(x) x(end-7:end), cellstr(fileNames),...
        'UniformOutput', false);
    [~, idx] = sort(shortNames);
    fileNames = fileNames(idx);

    numFiles = numel(fileNames);
    fprintf('Found %d files in %s\n', numFiles, fPath);

    % Read in each file
    data = cell(numFiles, 1);
    for i = 1:numFiles
        data{i} = dlmread(fullfile(fPath, fileNames(i)));
    end

