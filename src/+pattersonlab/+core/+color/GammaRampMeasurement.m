classdef GammaRampMeasurement < handle
% GAMMARAMPMEASUREMENT
%
% Description:
%   Quantifies a light source's input-output function as measured by an
%   Ocean Optics spectrometer. Assumes one of the measurements is at 0 to
%   provide a noise floor.
%
% Constructor:
%   obj = GammaRampMeasurement("420nm", "C:\420nm", 5:-0.1:0, 2, "20220422")
%
% Inputs:
%   lightSource             string
%       Light source's name (e.g., "420 nm" or "Blue")
%   filePath                string, must be a folder
%       The folder where the measurements are located
%   values                  vector
%       The light source's inputs for each measurement (e.g., voltage to an
%       LED or 8-bit value to a projector). These must match the order in
%       which the measurements were collected (by default, the first of the
%       numbers automatically added after "AbsoluteIrradiance__". See Token
%   beamDiameter            double
%       Diameter of beam on sensor in mm
%   calibrationDate         string
%       Date calibration was performed (I use YYYYmmDD but don't enforce
%       conversion to datetime since we mainly use it for file naming)
% Optional key/value inputs:
%   Token                   string
%       A token used to identify the desired measurement .txt files. Use to
%       omit unrelated .txt files (default = "_AbsoluteIrradiance_")
%   ND                      double
%       Value of neutral density filter(s) in place when measurement was
%       made (default is 0)
%   MinWL                   double  (default = N/A)
%       Minimum wavelength expected for spectra (values for wavelengths
%       below this will be set to zero to minimize impact of noise)
%   MaxWL                   double  (default = N/A)
%       Maximum wavelength expected for spectra (values for wavelengths
%       above this will be set to zero to minimize impact of noise)
%
% Use:
%   obj = GammaRampMeasurement("420nm", "C:\420nm", 5:-0.1:0, 2, "20220422",...
%       "Token", "420nm_Abs", "MinWL", 400, "MaxWL", 450);
%   % Plot the normalized spectra, then write to file
%   obj.plotSpectra(); obj.writeSpectra("C:\..\MySpectraFolder");
%   % See the input-output lookup table (optionally upsample values)
%   obj.plotLUT(0:0.05:5);
%   % Save the lookup table
%   obj.writeLUT("C:\..\MySpectraFolder", )
%
% History:
%   7Apr2024 - SSP - Adapted from Tyler's script
%   12Dec2024 - SSP - Debugged some table outputs
%   03Apr2025 - SSP - Added NDF, fit functions, and omitted values
%   09Apr2025 - SSP - Added JSON export and import
% --------------------------------------------------------------------------

    properties (SetAccess = private)
        lightSource                 (1,1)       string
        values                      (1,:)       {mustBeNumeric}
        Measurements
        calibrationDate
        beamDiameter                (1,1)       double       % mm
        NDF                         (1,1)       double = 0
        omittedValues                           double = []
        importProps                 (1,1)       struct
        lutFitType                  (1,1)       string = "none"
    end

    properties (Dependent)
        numMeasurements             (1,1)       {mustBeInteger}
        intensities                             double
        spectra                                 double
        wavelengths                             double
        powers                                  double
        maxPower                    (1,1)       double
        peakWavelength              (1,1)       double
    end

    properties (Hidden, Constant)
        DEFAULT_OUTPUT_FOLDER = fullfile(pattersonlab.core.color.util.getMainFolder(), "output");
    end

    methods
        function obj = GammaRampMeasurement(lightSource, filePath, ...
                values, beamDiameter, calibrationDate, opts)
            arguments
                lightSource     (1,1)   string
                filePath        (1,1)   string
                values          (1,:)   double
                beamDiameter    (1,1)   double      {mustBePositive}
                calibrationDate (1,1)   string
                opts.Token      (1,1)   string      = "_AbsoluteIrradiance_"
                opts.NDF        (1,1)   double      {mustBeNonnegative} = 0
                opts.MinWL              double      {mustBeScalarOrEmpty} = []
                opts.MaxWL              double      {mustBeScalarOrEmpty} = []
                opts.OmittedValues                  = []
                opts.FitType    (1,1)   string      = "none"
                opts.ImportData (1,1)   logical     = true
            end

            obj.lightSource = lightSource;
            obj.values = values;
            obj.NDF = opts.NDF;
            obj.beamDiameter = beamDiameter;
            obj.calibrationDate = calibrationDate;
            obj.setOmittedValues(opts.OmittedValues);

            obj.importProps = struct(...
                'Folder', filePath, 'Token', opts.Token,...
                'MinWL', opts.MinWL, 'MaxWL', opts.MaxWL);

            if opts.ImportData
                obj.importMeasurements();
            end
        end

        function addMeasurement(obj, measurement)
            % ADDMEASUREMENT
            %
            % Description:
            %   Append a new measurement (for loading from a file)
            % --------------------------------------------------------------
            arguments
                obj
                measurement         pattersonlab.core.color.SpectralMeasurement
            end

            obj.Measurements = [obj.Measurements; measurement];
        end

        function setOmittedValues(obj, badValues)
            if isempty(badValues)
                obj.omittedValues = [];
                return
            end

            mustBeMember(badValues, obj.values);
            obj.omittedValues = sort(badValues);
        end

        function setNDF(obj, ndf)
            if isempty(ndf)
                ndf = 0;
            end
            obj.NDF = ndf;
            obj.refresh();
        end

        function setMinWL(obj, minWL)
            obj.importProps.MinWL = minWL;
            obj.refresh();
        end

        function setMaxWL(obj, maxWL)
            obj.importProps.MaxWL = maxWL;
            obj.refresh();
        end

        function setFitType(obj, fitType)
            arguments
                obj
                fitType (1,1) string {mustBeMember(fitType, ["none", "linear", "poly8"])}
            end

            obj.lutFitType = fitType;
        end

        function refresh(obj)
            obj.importMeasurements();
        end
    end

    % Dependent set/get methods
    methods
        function out = get.numMeasurements(obj)
            out = numel(obj.Measurements);
        end

        function out = get.intensities(obj)
            out = arrayfun(@(x) x.voltage, obj.Measurements);
            out = round(out, 10);
        end

        function out = get.wavelengths(obj)
            out = obj.Measurements(end).wavelengths;
        end

        function out = get.powers(obj)
            out = arrayfun(@(x) x.power, obj.Measurements);
            if obj.Measurements(1).voltage == 0
                out(1) = 0;
            end
        end

        function out = get.maxPower(obj)
            out = max(obj.powers);
        end

        function out = get.spectra(obj)
            out = obj.getNormalizedSpectra();
        end

        function out = get.peakWavelength(obj)
            out = obj.Measurements(end).peakWavelength();
        end
    end

    % Core methods
    methods
        function out = getNormalizedSpectra(obj, targetValue)
            if nargin < 2
                targetValue = max(obj.intensities);
            end

            idx = find(obj.intensities == targetValue);
            out = obj.Measurements(idx).spectra;
            out = out / max(out);
        end

        function T = getLUT(obj, targetValues, fitType)
            if nargin < 3
                fitType = obj.lutFitType;
            else
                obj.lutFitType = fitType;
            end

            if any(~ismember(targetValues, obj.intensities))
                % Interpolate to the provided target values if needed
                pwr = obj.interpolate(targetValues);
            else % Use the default values provided
                pwr = obj.powers;
            end

            [~, idx] = sort(targetValues);

            if isrow(targetValues)
                targetValues = targetValues';
            end
            if isrow(pwr)
                pwr = pwr';
            end

            T = table(targetValues(idx), pwr(idx),...
                'VariableNames', {'Input', 'Power'});
            T = obj.fitLut(T, fitType, targetValues);
        end
    end

    % File output methods
    methods
        function fName = getFileName_LUT(obj)
            fName = sprintf('LUT_%s_%s_%sndf.txt',...
                obj.lightSource, obj.calibrationDate, num2str(10*obj.NDF));
        end

        function fName = getFileName_Spectra(obj)
            fName = sprintf('%s_%s_%sndf.txt',...
                    obj.lightSource, obj.calibrationDate, num2str(10*obj.NDF));
        end

        function writeLUT(obj, savePath, targetValues, fitType)
            % Saves the lookup table in a file matching Qiang's AOSLO
            % software specifications (e.g., 0.1 V increments)
            arguments
                obj
                savePath        (1,1)   string  {mustBeFolder} = obj.DEFAULT_OUTPUT_FOLDER
                targetValues    (1,:)   double  = 0:0.1:5
                fitType         (1,1)   string  = obj.lutFitType
            end

            fName = obj.getFileName_LUT();
            fName = fullfile(savePath, fName);

            T = obj.getLUT(targetValues, fitType);
            % Variable names wanted by Qiang's AOSLO code...
            T.Properties.VariableNames = {'VOLTAGE', 'POWER'};
            writetable(T, fName, 'Delimiter', '\t')
            fprintf('Saved LUT as %s\n', fName);
        end

        function writeSpectra(obj, savePath, targetValue, fName)
            % WRITESPECTRA
            %
            % Description:
            %   Write the spectra to a .txt file (for use in deprecated
            %   calibration code, no longer necessary)
            % -------------------------------------------------------------
            arguments
                obj
                savePath        (1,1)   string  {mustBeFolder} = obj.DEFAULT_OUTPUT_FOLDER
                targetValue     (1,1)   double  = max(obj.intensities)
                fName           (1,1)   string  = ""
            end

            data = obj.getNormalizedSpectra(targetValue);
            if fName == ""
                fName = sprintf('%s_%s_%sndf.txt',...
                    obj.lightSource, obj.calibrationDate, num2str(10*obj.NDF));
            end

            fName = fullfile(savePath, fName);
            T = table(obj.wavelengths, data,...
                'VariableNames', {'Lambda', 'Normalized'});
            writetable(T, fName, 'Delimiter', '\t')
            fprintf('Saved normalized spectra as %s\n', fName);
        end

        function exportToJson(obj, fName)
            fName = convertCharsToStrings(fName);
            if ~endsWith(fName, ".json")
                fName = fName + ".json";
            end
            if isempty(fileparts(fName))
                outputDir = fullfile(pattersonlab.core.color.util.getMainFolder(), "output");
                fName = fullfile(outputDir, fName);
            end

            S = struct(obj);

            writestruct(S, fName);
            fprintf('Wrote file: %s\n', fName);
        end
    end

    % Plotting methods
    methods
        function plotSpectra(obj, opts)
            arguments
                obj
                opts.Parent                         = []
                opts.Area       (1,1)   logical     = true
                opts.Color                          = [0.5 0.5 1]
                opts.ShowCutoffs (1,1)  logical     = false
            end

            if isempty(opts.Parent)
                axHandle = axes('Parent', figure('Name', 'NormSpectra'));
                hold on;
                xlim([380, 720]); xlabel('Wavelength (nm)');
                ylim([0 1]); ylabel('Normalized');
                title(obj.lightSource);
                grid on;
            else
                axHandle = opts.Parent;
            end

            if opts.Area
                area(obj.wavelengths, obj.spectra,...
                    "FaceColor", opts.Color, "EdgeColor", [0.1 0.1 0.1],...
                    "LineWidth", 1, "FaceAlpha", 0.5, "Parent", axHandle);
            else
                plot(axHandle, obj.wavelengths, obj.spectra,...
                    'Color', opts.Color, 'LineWidth', 1);
            end
            if opts.ShowCutoffs
                plot([obj.importProps.MinWL obj.importProps.MinWL], [0 1],...
                    'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0.25 0.25]);
                plot([obj.importProps.MaxWL obj.importProps.MaxWL], [0 1],...
                    'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0.25 0.25]);
            end
        end

        function plotLUT(obj, targetValues, opts)
            arguments
                obj
                targetValues        (1,:)   double  = obj.intensities
                opts.Parent                         = []
                opts.Color                          = [0 0 1]
                opts.FitType        (1,1)   string  = obj.lutFitType
            end

            lut = obj.getLUT(targetValues, opts.FitType);

            if isempty(opts.Parent)
                opts.Parent = axes('Parent', figure('Name', sprintf('%s LUT', obj.lightSource)));
            end

            hold(opts.Parent, 'on');
            plot(lut.Input, lut.Power, '-o',...
                'MarkerEdgecolor', opts.Color,...
                'MarkerFaceColor', lighten(opts.Color, 0.5));
            xlabel('Voltage (V)'); ylabel('Power (uW)');
        end

        function plotBoth(obj, targetValues)
            if nargin < 2
                targetValues = obj.intensities;
            end

            figure();
            subplot(2, 1, 1); hold on;
            obj.plotSpectra("Parent", gca);
            subplot(2, 1, 2);
            obj.plotLUT(targetValues, "Parent", gca);
        end

        function compareNormalizedSpectra(obj, plotValues)
            obj.plotSpectra('Area', false);

            allValues = obj.intensities;

            for i = 1:numel(plotValues)
                if plotValues(i) == max(obj.intensities)
                    continue
                end
                idx = find(allValues == plotValues(i));
                if isempty(idx)
                    error('Value %s not found', num2str(plotValues(i)));
                end
                plot(obj.Measurements(idx).wavelengths, obj.getNormalizedSpectra(plotValues(i)));
            end
        end
    end

    methods (Access = private)
        function importMeasurements(obj)
            [data, fileNames] = pattersonlab.core.color.io.loadSpectralMeasurementFiles(...
                obj.importProps.Folder, obj.importProps.Token);
            if numel(obj.values) ~= numel(fileNames)
                error('importMeasurements:InvalidFileCount',...
                    '%u values were provided, but %u files were found',...
                    numel(obj.values), numel(fileNames));
            end
            if obj.NDF == 0
                ndfString = [];
            else
                ndfString = arrayfun(@(x) sprintf("NE%sA-A.txt", int2fixedwidthstr(10*x, 2)), obj.NDF,...
                    "UniformOutput", true);
            end

            obj.Measurements = [];
            for i = 1:numel(fileNames)
                newMeasurement = pattersonlab.core.color.SpectralMeasurement(data{i}, obj.values(i),...
                    obj.beamDiameter, obj.importProps.MinWL, obj.importProps.MaxWL,...
                    ndfString);
                obj.Measurements = [obj.Measurements, newMeasurement];
            end
            obj.Measurements = obj.Measurements';

            [sortedValues, idx] = sort(obj.values);
            obj.Measurements = obj.Measurements(idx);

            if sortedValues(1) == 0
                arrayfun(@(x) setBackground(x, obj.Measurements(1).spectra0), obj.Measurements(2:end));
            end

        end

        function out = interpolate(obj, targetValues)
            if any(targetValues < min(obj.intensities)) || any(targetValues > max(obj.intensities))
                error('Values out of range');
            end

            validIntensities = obj.intensities(~ismember(obj.intensities, obj.omittedValues));
            validPowers = obj.powers(~ismember(obj.intensities, obj.omittedValues));

            out = interp1(validIntensities, validPowers, targetValues);
        end
    end

    % MATLAB built-in methods
    methods
        function S = struct(obj)
            S = struct(...
                'lightSource', obj.lightSource,...
                'values', obj.values,...
                'beamDiameter', obj.beamDiameter,...
                'calibrationDate', obj.calibrationDate,...
                'NDF', obj.NDF,...
                'omittedValues', obj.omittedValues,...
                'LutFitType', obj.lutFitType,...
                'importProps', obj.importProps);
            S.Measurements = arrayfun(@(x) x.struct(), obj.Measurements);
        end
    end

    methods (Static)
        function obj = initFromJSON(S)
            % Create from a struct loaded from a JSON file
            if ~isstruct(S)
                S = readstruct(S);
            end

            obj = pattersonlab.core.color.GammaRampMeasurement(...
                S.lightSource, S.importProps.Folder,...
                S.values, S.beamDiameter, S.calibrationDate,...
                'Token', S.importProps.Token, 'NDF', S.NDF,...
                'MinWL', S.importProps.MinWL, 'MaxWL', S.importProps.MaxWL,...
                "ImportData", false, 'FitType', S.LutFitType);

            obj.setOmittedValues(S.omittedValues);

            for i = 1:numel(S.Measurements)
                m = pattersonlab.core.color.SpectralMeasurement.initFromJSON(S.Measurements(i));
                obj.addMeasurement(m);
            end
        end
    end

    methods (Static, Access = private)
        function lutFit = fitLut(lut, fitType, targetValues)

            targetValues = targetValues(:);

            switch fitType
                case 'poly8'
                    [fitFcn, gof] = fit(lut.Input, lut.Power,...
                        fittype('poly8'), 'Normalize', 'on');
                case 'linear'
                    ft = fittype('a*x+b', 'Independent', 'x', 'Dependent', 'y');
                    opts = fitoptions('Method', 'NonlinearLeastSquares',...
                        'StartPoint', [max(lut.Power)/max(lut.Input) 0]);

                    [fitFcn, gof] = fit(lut.Input, lut.Power, ft, opts);
                case 'none'
                    lutFit = lut;
                    return
            end

            lutFit = table(targetValues, fitFcn(targetValues),...
                'VariableNames', {'Input', 'Power'});
            lutFit = sortrows(lutFit, "Input");
            lutFit.Power(1) = lut.Power(1);
            lutFit.Power(end) = lut.Power(end);
            lutFit.Power(lutFit.Power < 0) = 0;

            fprintf('Fit %s r2 = %.3f\n', fitType, gof.rsquare);
        end
    end
end