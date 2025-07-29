function X = makeLEDStimulusFile(fName, X, lutFiles)
% MAKELEDSTIMULUSFILE
%
% Syntax:
%   X = makeLEDStimulusFile(fName, X);
%
% Inputs:
%   fName                       char
%       File name (ending with .txt). If a full file path isn't
%       specified, the file will save in a folder called "stimuli"
%   X                           [N x 3]
%       LED powers at each 2 ms time point
%   lutFiles                    string (1 x 3)
%       File names for each LED's LUT
%
% History:
%   08Dec2021 - SSP
%   03Nov2022 - SSP - Added calibration obj for LUT file names
%   19Apr2025 - SSP - More generic input for LUT file names
% --------------------------------------------------------------------------

    % Input validation
    fName = convertCharsToStrings(fName);
    if ~endsWith(fName, '.txt')
        fName = fName + ".txt";
    end
    if fileparts(fName) == ""
        outputDir = fullfile(pattersonlab.core.color.util.getMainFolder(), "output");
        fName = fullfile(outputDir, fName);
    end

    if size(X, 1) <= 3
        X = X';
    end
    assert(size(X, 2) == 3, 'X must be an [N x 3] matrix!');

    if isa(lutFiles, 'pattersonlab.core.color.LedCalibration')
        obj = lutFiles;
        lutFiles = [...
            string(obj.RED.getFileName_LUT()),...
            string(obj.GREEN.getFileName_LUT()),...
            string(obj.BLUE.getFileName_LUT())];
    end

    % Open and discard existing contents
    fid = fopen(fName, 'w');
    if fid == -1
        error('File %u could not be opened!', fName);
    end

    % Write the stimulus metadata
    fprintf(fid, '[header]\r\n');
    fprintf(fid, 'functionality    = 1\r\n');

    % Add the LUT file names (note: this is an absolute path on imaging PC)
    fprintf(fid,...
        'lut1 		= C:\\Users\\slo\\Stimuli\\LUTs\\%s\r\n', lutFiles(1));
    fprintf(fid,...
        'lut2 		= C:\\Users\\slo\\Stimuli\\LUTs\\%s\r\n', lutFiles(2));
    fprintf(fid,...
        'lut3 		= C:\\Users\\slo\\Stimuli\\LUTs\\%s\r\n', lutFiles(3));

    fprintf(fid, 'interval_value	= 2\r\n');
    fprintf(fid, 'interval_unit	= ms\r\n');
    fprintf(fid, ['data_len	=', num2str(size(X, 1)) '\r\n']);

    % Write the stimulus data
    fprintf(fid, '\r\n');
    fprintf(fid, '[data]\r\n');
    for i = 1:size(X, 1)
        fprintf(fid, '%u=%.4f,%.4f,%.4f\r\n', i, X(i, 1), X(i, 2), X(i, 3));
    end

    fclose(fid);

    % Report to the command line
    fprintf('\tCompleted %s\\%s\n', fName);