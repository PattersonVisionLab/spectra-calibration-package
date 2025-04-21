classdef LedCalibration < handle
% LEDCALIBRATION
%
%
% Key methods
%   writeSpectra(obj, lutDir)
%   saveToJSON(obj, savePath)
% --------------------------------------------------------------------------

    properties (SetAccess = private)
        calibrationDate

        RED                 % LED #1
        GREEN               % LED #2
        BLUE                % LED #3

        NDF

        wavelengths         % Wavelength array (nm)

        lmsLambdaMax        % Peak wavelengths of the L, M and S cones
        cones               % L, M, S-cones spectra, ready for calcs
        receptors           % Cones and luminance

        primaries           % LED spectra preprocessed, ready for calcs
        ledPowers           % Max powers for each LED (uW)
        ledMeans            % Mean values for each LED (uW)

        ndfAttenuations

        meanChromaticity    % CIE 1931 coordinates of white point
        stimPowers          % Table of LED powers for each stimulus
        stimContrasts       % Cone contrasts of each stimulus
    end

    % Properties used in calculations, usually not needed afterwards
    properties (Hidden, SetAccess = private)
        transmittance
        absorptance
        lum

        relativePowers
        spectraSplined
        rawSpectra
        ambientSpectra

        ledDir
        ledFiles
        lutFiles

        xyz
        primaryHeadRoom
        maxPowerDiff
        includeMacularPigment

        autoPlot
        verbose
    end

    properties (Dependent)
        numWavelengths
        meanPower           % Total power of LED means (uW)
        ledBkgd             % LED powers at 50%
    end

    properties (Hidden, Constant)
        % LEDs and cones
        NUM_LEDS = 3;
        CONES = {'L','M','S'};
        DEFAULT_LAMBDA_MAX = [660 530 420];
        LED_PLOT_COLORS = [1 0.25 0.25; 0 0.8 0.3; 0.2 0.3 0.9];

        % Misc parameters
        DEFAULT_HEADROOM  = 0.002;
        DARK_VALUES = [0 0 0];
    end



    methods
        function obj = LedCalibration(calibrationData, varargin)
            % LEDCALIBRATION
            %
            % ---------------------------------------------------------

            obj.RED = calibrationData(1);
            obj.GREEN = calibrationData(2);
            obj.BLUE = calibrationData(3);
            obj.ledPowers = [obj.RED.maxPower, obj.GREEN.maxPower, obj.BLUE.maxPower];

            ip = inputParser();
            ip.CaseSensitive = false;
            addParameter(ip, 'Bkgd', [], @isnumeric);
            addParameter(ip, 'MacularPigment', true, @islogical);
            addParameter(ip, 'LMSLambdaMax', obj.DEFAULT_LAMBDA_MAX, @isnumeric);
            addParameter(ip, 'Verbose', true, @islogical);
            addParameter(ip, 'Plot', false, @islogical);
            parse(ip, varargin{:});

            obj.ledMeans = ip.Results.Bkgd;

            obj.verbose = ip.Results.Verbose;
            obj.lmsLambdaMax = ip.Results.LMSLambdaMax;
            obj.includeMacularPigment = ip.Results.MacularPigment;
            obj.autoPlot = ip.Results.Plot;

            obj.NDF = unique([obj.RED.NDF, obj.GREEN.NDF, obj.BLUE.NDF]);
            obj.calibrationDate = obj.RED.calibrationDate;

            % Input checking
            if size(obj.ledMeans, 1) == 1
                obj.ledMeans = obj.ledMeans';
            end

            obj.initialize();
        end
    end

    % Dependent set/get methods
    methods
        function meanPower = get.meanPower(obj)
            meanPower = sum(obj.ledMeans .* obj.ledPowers');
        end

        function numWavelengths = get.numWavelengths(obj)
            numWavelengths = numel(obj.wavelengths);
        end

        function value = get.ledBkgd(obj)
            T = obj.stimPowers('Lum');
            value = T.Bkgd;
        end
    end

    % Plotting methods
    methods
        function ax = plotCIE(obj)
            % PLOTCIE
            %
            % Syntax:
            %   ax = obj.plotCIE()
            % ---------------------------------------------------------
            ax = axes('Parent', figure());
            plotChromaticity();
            hold on;
            plot(obj.meanChromaticity(1), obj.meanChromaticity(2),...
                'xk', 'MarkerSize', 10, 'LineWidth', 1, 'Tag', 'Bkgd');
            xy1 = [obj.getCieCoords([obj.ledMeans(1) 0 0])];
            xy2 = [obj.getCieCoords([0 obj.ledMeans(2) 0])];
            xy3 = [obj.getCieCoords([0 0 obj.ledMeans(3)])];
            plot([xy1(1) xy2(1) xy3(1) xy1(1)], [xy1(2) xy2(2) xy3(2) xy1(2)],...
                'k--', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'Tag', 'LEDs');
            figPos(gcf, 0.7, 0.7);
        end
    end

    % Analysis methods
    methods
        function initialize(obj)

            % Receptors matrix for L, M, S and luminance
            obj.initCones();
            obj.initLuminance();
            obj.receptors = [obj.cones; obj.lum];

            % CIE 1931 XYZ
            obj.initXYZ();

            % LED calibrations
            obj.initPrimaries();

            % Ambient light
            obj.initAmbientLight();

            % Background
            if isempty(obj.ledMeans)
                obj.ledMeans = obj.fitBackground()';
            end
            obj.initBackground();

            % Psychtoolbox parameters that we aren't currently using
            obj.initMisc();

            % Where stimulus information will be saved
            obj.stimContrasts = containers.Map();
            obj.stimPowers = containers.Map();

            % Calculate initial cone isolating stimuli
            obj.runConeIsolation();
        end

        function runConeIsolation(obj)
            % RUNCONEISOLATION
            %
            % Syntax:
            %   obj.runConeIsolation()
            % ---------------------------------------------------------
            disp('Hello ------------------------')
            import patterson.core.color.lib.*
            names = {'L', 'M', 'S', 'LM', 'Isolum', 'Lum'};
            targets = {1, 2, 3, [1 2], [1 2], 1:4, 1:4};
            contrasts = {[], [], [], [], [40 40], 0.95+ones(1, 4)};
            nReceptors = [3, 3, 3, 3, 4, 4, 4];

            % SST parameters not used for these calculations
            receptorsToIgnore = [];
            receptorsToMinimize = [];
            primariesToPin = [];

            for i = 1:numel(names)
                % Receptor isolating.  Specify which receptor you want to
                % isolate and to ignore luminance. Here what works best is
                % to get rid of the luminance row in the receptors matrix
                receptorsToTarget = targets{i};
                desiredContrast = contrasts{i};

                % Careful examaination of the arguments will reveal that
                % the initialGuess for the primaries is set to the
                % background value for the primaries, so that the
                % constraints are all met at the start of the search.  The
                % optimization routine is much happier when it is done this
                % way -- bad things happen if you start with a guess that
                % violates constraints.
                positiveModulationPrimary = ReceptorIsolate(...
                    obj.receptors(1:nReceptors(i),:),...
                    receptorsToTarget, receptorsToIgnore, receptorsToMinimize,...
                    obj.primaries, obj.ledMeans, obj.ledMeans, primariesToPin,...
                    obj.primaryHeadRoom, obj.maxPowerDiff, desiredContrast,...
                    obj.ambientSpectra);

                % For a sinusoidal modulation, you'd multipy deltaPrimary
                % by the sinusoidal value (ranging between -1 and 1) and
                % add this to the background primary
                deltaPrimary = positiveModulationPrimary - obj.ledMeans;
                negativeModulationPrimary = obj.ledMeans - deltaPrimary;

                % Luminance is included here. That just lets us see
                % luminance contrast even if we didn't specify it when we
                % computed.  Note that L and M cone isolating modulations
                % also produce luminance contrast.
                backgroundReceptors = obj.receptors * (obj.primaries * obj.ledMeans + obj.ambientSpectra);
                deltaReceptors = obj.receptors * obj.primaries * (positiveModulationPrimary - obj.ledMeans);
                contrastReceptors = deltaReceptors ./ backgroundReceptors;

                % Print results to command line
                if obj.verbose
                    fprintf('%s contrast: %s\n', names{i}, ...
                        num2str(round(contrastReceptors(targets{i}),3)'));
                end

                obj.stimPowers(names{i}) = table(...
                    obj.ledMeans .* obj.ledPowers',...
                    positiveModulationPrimary .* obj.ledPowers',...
                    negativeModulationPrimary .* obj.ledPowers',...
                    deltaPrimary .* obj.ledPowers',...
                    'RowNames', {'R', 'G', 'B'},...
                    'VariableNames', {'Bkgd', 'Up', 'Down', 'dP'});

                obj.stimContrasts(names{i}) = contrastReceptors;
            end

            if obj.verbose
                fprintf('White point power at 50%% = %.3f uW\n', obj.meanPower);
            end
        end

        function ledScalars = fitBackground(obj)
            % FITBACKGROUND
            %
            % Description:
            %   Identifies LED weights that achieve 0.33 0.33 CIE
            %
            % Syntax:
            %   ledScalars = obj.fitBackground()
            % ---------------------------------------------------------
            objFcn = @(powers) norm(obj.getCieCoords(powers) - [0.33 0.33]);

            options = optimoptions('fmincon', ...
                'Display', 'final', 'Algorithm', 'sqp');

            lb = [0 0 0];
            ub = obj.ledPowers;
            initialPowers = 0.5 * obj.ledPowers;

            scaledPowers = fmincon(objFcn, initialPowers, [], [], [], [],...
                lb, ub, [], options);

            ledScalars = scaledPowers ./ obj.ledPowers;
            cieCoords = obj.getCieCoords(scaledPowers);

            if min(ledScalars) > 0.5
                ledScalars = ledScalars/median(ledScalars)/2;
            end

            fprintf('White Point fit: %.3f %.3f for values %.3f %.3f %.3f\n',...
                cieCoords, ledScalars);
        end

        function xy = getCieCoords(obj, ledPowers)
            import patterson.core.color.lib.*
            scalars = ledPowers ./ obj.ledPowers;
            backgroundXYZ = obj.xyz * obj.primaries * scalars';
            backgroundxyY = XYZToxyY(backgroundXYZ);
            xy = backgroundxyY(1:2)';
        end

        function stim = calcStimulus(obj, whichStim, baseStim)
            if isa(whichStim, 'sara.SpectralTypes')
                whichStim = whichStim.getAbbrev();
            end

            T = obj.stimPowers(whichStim);
            dPower = T.dP';
            bkgdPower = T.Bkgd';

            try
                stim = (dPower .* ((1/baseStim(1)) * ...
                    (baseStim-baseStim(1))) + bkgdPower);
            catch
                stim = zeros(3, numel(baseStim));
                for i = 1:3
                    stim(i,:) = (dPower(i) .* ((1/baseStim(1)) * ...
                        (baseStim-baseStim(1))) + bkgdPower(i));
                end
            end
        end
    end


    % Initialization methods
    methods (Access = private)
        function initCones(obj)
            import pattersonlab.core.color.lib.*;
            % INITCONES
            obj.wavelengths = (380:750)';
            S = WlsToS(obj.wavelengths);
            lcone = StandardTemplate(obj.lmsLambdaMax(1), 350:700);
            mcone = StandardTemplate(obj.lmsLambdaMax(2), 350:700);
            scone = StandardTemplate(obj.lmsLambdaMax(3), 350:700);

            obj.cones = SplineCmf((350:700)', [lcone; mcone; scone], S);
            obj.transmittance = importdata('rhesus_transmittance.txt', '\t');
            if obj.includeMacularPigment
                obj.absorptance = macpigment(obj.wavelengths)';
            else
                obj.absorptance = zeros(length(obj.wavelengths), 1);
            end

            % Convert transmittance into a percentage
            obj.transmittance = interp1(...
                obj.transmittance(:,1)', obj.transmittance(:,2)',...
                obj.wavelengths, 'spline', 'extrap');
            obj.transmittance = obj.transmittance ./ 100;
            obj.transmittance(obj.transmittance < 0) = 0;
            obj.absorptance = 1 - interp1(obj.wavelengths', obj.absorptance',...
                obj.wavelengths, 'spline', 'extrap');

            % Apply to cone spectral sensitivities
            obj.cones = obj.cones .* obj.transmittance' .* obj.absorptance';
            obj.cones = obj.cones ./ max(obj.cones, [], 2);
        end

        function initXYZ(obj)
            % CIE 1931 XYZ
            import pattersonlab.core.color.lib.*;
            load('T_xyz1931.mat'); %#ok<LOAD>
            obj.xyz = SplineCmf(S_xyz1931, T_xyz1931, WlsToS(obj.wavelengths));
        end

        function initLuminance(obj)
            % Specify luminance as weighted sum of L and M cones.  We want
            % this to be a weighted sum of the L and M cones, as that
            % matches how we typically think about post-receptoral
            % channels.  To enforce this, we start with a luminance
            % spectral sensitivity and fit it with a weigthed sum of the L
            % and M cones that we're using.
            vlambda = obj.cones(1,:) + obj.cones(2,:);
            % Defines luminance as 1:1 M:L for macaque
            lumTabulated = vlambda ./ max(vlambda);
            lumFactorsLM = (obj.cones(1:2, :)' \ lumTabulated');
            obj.lum = (obj.cones(1:2, :)' * lumFactorsLM)';
        end

        function initAmbientLight(obj)
            % The ambient light is not the background for the modulations.
            % Rather it is the light delivered to the eyes when the
            % settings on the primaries are all zero.  With LED's, this may
            % well be zero.  With typical computer displays there is often
            % sone light even when all the DAC values are 0, which is why
            % we allow specification of something here.
            %
            % This is specified as a spectrum, as in our lab we typically
            % measure it directly and sometimes it is not a combination of
            % the primarie (for example, if there is ambient in the room
            % that that comes from some other source).
            %
            % Currently we assume it is dark when all primaries are set to
            % 0, and this is probably fine as an approximation.
            obj.ambientSpectra = zeros(size(obj.wavelengths));
        end

        function initPrimaries(obj)
            obj.primaries = zeros(numel(obj.wavelengths), numel(obj.ledFiles));

            obj.rawSpectra = cell(1, 3);
            obj.relativePowers = zeros(1, 3);
            obj.spectraSplined = cell(1,3);
            ledSpectraRaw = cell(1, 3);

            ledSpectraRaw{1} = [obj.RED.wavelengths, obj.RED.getNormalizedSpectra()];
            ledSpectraRaw{2} = [obj.GREEN.wavelengths, obj.GREEN.getNormalizedSpectra()];
            ledSpectraRaw{3} = [obj.BLUE.wavelengths, obj.BLUE.getNormalizedSpectra()];

            for i = 1:obj.NUM_LEDS
                obj.rawSpectra{i} = ledSpectraRaw{i};
                obj.spectraSplined{i} = interp1(...
                    ledSpectraRaw{i}(:,1), ledSpectraRaw{i}(:,2),...
                    obj.wavelengths, 'linear', 0);

                % Get constant by integrating relative spectra over
                % wavelength. Psychtoolbox, which underlies this code,
                % likes to "think" of power in units of power/wl-band
                % rather than power/nm and we follow that convention here.
                % If we are operating on 1 nm wavelength sampling, as we
                % were when this code was first written, the two
                % conventions collapse to the same thing.  Be a little
                % careful about this if you change the way you bring in
                % your primary spectra.
                obj.relativePowers(i) = sum(obj.spectraSplined{i});
            end

            % Separate so it's easy to trial different LED powers
            obj.scaleByPowerMeasurement();
        end

        function scaleByPowerMeasurement(obj, ledPowers)
            if nargin == 2
                obj.ledPowers = ledPowers;
                disp('Set new ledPowers');
            end

            obj.primaries = zeros(numel(obj.wavelengths), obj.NUM_LEDS);
            for i = 1:obj.NUM_LEDS
                obj.primaries(:,i) = obj.ledPowers(i) * obj.spectraSplined{i} / obj.relativePowers(i);
            end
        end

        function initBackground(obj)
            % Calculate the chromaticity (should be 0.33 0.33)
            import pattersonlab.core.color.lib.*;

            backgroundXYZ = obj.xyz * obj.primaries * obj.ledMeans;
            obj.meanChromaticity = XYZToxyY(backgroundXYZ);
            fprintf('Background 1931 xy chromaticity: %0.3f, %0.3f\n\n',...
                obj.meanChromaticity(1:2));
        end

        function initMisc(obj)
            % This really doesn't do anything, leave it as a big number.
            obj.maxPowerDiff = 10000;

            % This parameter says how close to the edge of the gamut we
            % allow a modulation to get.  So if we specify 0.02, the
            % primary values can range between 0.02 and 0.98. I like to
            % leave a little headroom, becuase devices can get a little
            % wonky right at their extrema, and sometimes the gamut change
            % over time and I like to avoid having to change stimuli if the
            % gamut gets a little smaller.
            if isempty(obj.primaryHeadRoom)
                obj.primaryHeadRoom = obj.DEFAULT_HEADROOM;
            end
        end
    end

    % Input/output methods
    methods
        function writeSpectra(obj, lutDir)
            arguments
                obj
                lutDir      {mustBeFolder}
            end

            obj.RED.writeSpectra(lutDir);
            obj.GREEN.writeSpectra(lutDir);
            obj.BLUE.writeSpectra(lutDir);
        end

        function exportToJSON(obj, fPath)
            if nargin < 2
                fPath = fullfile(pattersonlab.core.color.util.getMainFolder(), "output");
            end

            fName = sprintf("LedCalibration_%s_%undf.json",...
                obj.calibrationDate, 10*sum(obj.NDF));
            fName = fullfile(fPath, fName);

            S = struct(obj);

            writestruct(S, fName);
            fprintf('Wrote file: %s\n', fName);
        end
    end

    % MATLAB built-in methods
    methods
        function S = struct(obj)
            S = struct(...
                'RedLED', obj.RED.struct(),...
                'GreenLED', obj.GREEN.struct(),...
                'BlueLED', obj.BLUE.struct(),...
                'LmsLambdaMax', obj.lmsLambdaMax,...
                'MacularPigment', obj.includeMacularPigment,...
                'LedMaxPowers', obj.ledPowers,...
                'WhitePointWeights', obj.ledMeans,...
                'MeanChromaticity', obj.meanChromaticity);

            % Write the powers for each stimulus
            S.Stimuli = struct();
            S.Stimuli.Powers = struct();
            S.Stimuli.Contrasts = struct();
            k = obj.stimPowers.keys;

            % Background is always the same, only need it listed once
            T = obj.stimPowers(k{1});
            S.Stimuli.Powers.Background = round(T.Bkgd, 5)';

            for i = 1:numel(k)
                T = obj.stimPowers(k{i});
                S.Stimuli.Powers.(k{i}) = struct();
                S.Stimuli.Powers.(k{i}) = round(T.dP, 5)';
                S.Stimuli.Contrasts.(k{i}) = round(obj.stimContrasts(k{i}),5)';
            end

            S.Stimuli.Powers.Units = 'uW';
            S.Stimuli.Powers.Labels = {'R', 'G', 'B'};
            S.Stimuli.Powers.Description = 'The change in power from the background powers for a 100% contrast increment.';
            S.Stimuli.Contrasts.Units = 'Norm';
            S.Stimuli.Contrasts.Labels = {'L', 'M', 'S', 'Lum'};

            S.Files = struct(...
                'Spectra', [obj.RED.getFileName_Spectra(),...
                            obj.GREEN.getFileName_Spectra(),...
                            obj.BLUE.getFileName_Spectra()],...
                'LUT',     [obj.RED.getFileName_LUT(),...
                            obj.GREEN.getFileName_LUT(),...
                            obj.BLUE.getFileName_LUT()]);

            S.Spectra.Wavelengths = obj.wavelengths;
            S.Spectra.R = obj.RED.getNormalizedSpectra();
            S.Spectra.R(:,1) = [];
            S.Spectra.G = obj.GREEN.getNormalizedSpectra();
            S.Spectra.G(:,1) = [];
            S.Spectra.B = obj.BLUE.getNormalizedSpectra();
            S.Spectra.B(:,1) = [];

            % TODO: Hard-coded values
            rLUT = obj.RED.getLUT(5:-0.1:0);
            gLUT = obj.GREEN.getLUT(5:-0.1:0);
            bLUT = obj.BLUE.getLUT(5:-0.1:0);
            S.LUT = struct(...
                'Values', rLUT.Input,...
                'R', struct('Power', rLUT.Power,...
                    'FitType', obj.RED.lutFitType),...
                'G', struct('Power', gLUT.Power,...
                    'FitType', obj.GREEN.lutFitType),...
                'B', struct('Power', bLUT.Power,...
                    'FitType', obj.BLUE.lutFitType));
        end
    end



    methods (Static)
        function obj = initFromJSON(jsonFile)
            if ~isstruct(jsonFile)
                S = readstruct(jsonFile);
            else
                S = jsonFile;
            end

            import pattersonlab.core.color.*;

            red = GammaRampMeasurement.initFromJSON(S.RedLED);
            green = GammaRampMeasurement.initFromJSON(S.GreenLED);
            blue = GammaRampMeasurement.initFromJSON(S.BlueLED);

            obj = LedCalibration([red, green, blue],...
                'LmsLambdaMax', S.LmsLambdaMax,...
                'MacularPigment', S.MacularPigment,...
                'Bkgd', S.WhitePointWeights);
        end
    end
end