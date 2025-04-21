function figPos(fh, x, y)
% FIGPOS
%
% DESCRIPTION:
%   Change figure size by some factor
%
% SYNTAX:
%   figPos(fh, x, y)
%
% INPUTS:
%   fh    figure handle
%   x     factor to multiply width by
%   y     factor to multiply height by
%
% NOTES:
%   To keep x or y constant, input [] or 1
%
% History:
%   02Mar2017 - SSP
%   21Mar2017 - SSP - fixed screen position issue
%   14Oct2017 - SSP - removed unused pixel option, added misc checks
%   23Aug2020 - SSP - removed output, wasn't necessary, updated docs
% --------------------------------------------------------------------------

    if nargin < 3
        y = [];
    end

    pos = get(fh, 'Position');
    screenSize = get(0,'ScreenSize');

    if ~isempty(x) || x == 1
        pos(3) = pos(3) * x;
        if pos(3) > screenSize(3)
            % Keep figure size under screen size
            pos(3) = screenSize(3);
            pos(1) = screenSize(1);
            % Make sure figure isn't running off screen
        elseif pos(1) + pos(3) >= screenSize(3) - 50
            pos(1) = 50;
        end
    end

    if ~isempty(y) || y == 1
        pos(4) = pos(4) * y;
        if pos(4) > screenSize(4)
            pos(4) = screenSize(4);
            pos(2) = screenSize(2);
        elseif pos(2) + pos(4) > screenSize(4) - 50
            pos(2) = 50;
        end
    end

    set(fh, 'Position', pos);