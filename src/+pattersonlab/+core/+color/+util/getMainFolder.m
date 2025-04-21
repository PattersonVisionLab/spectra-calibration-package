function out = getMainFolder()
% GETMAINFOLDER
%
% Description:
%   Returns the path to the top-level folder for the repository.
%
% Syntax:
%   out = pattersonlab.core.color.io.getMainFolder()
%
% --------------------------------------------------------------------------

    out = fileparts(fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))));

