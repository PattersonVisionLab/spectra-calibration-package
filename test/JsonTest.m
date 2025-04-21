classdef JsonTest < matlab.unittest.TestCase

    properties
        TEST_PATH
        OUTPUT_PATH
        OBJECT
    end

    methods (TestClassSetup)
        function methodSetup(testCase)

            testCase.TEST_PATH = fullfile(pattersonlab.core.color.util.getMainFolder(), "test");
            testCase.OUTPUT_PATH = fullfile(testCase.TEST_PATH, "output");
            if ~exist(testCase.OUTPUT_PATH, "dir")
                mkdir(testCase.OUTPUT_PATH);
            end
            fprintf('Saving files to %s\n', testCase.OUTPUT_PATH);

            testCase.OBJECT = pattersonlab.core.color.io.loadMaxwellianViewSpectra(...
                fullfile(testCase.TEST_PATH, "data", "20250402_LEDSpectra"),...
                5:-0.2:0, "20250402", 3.9, 0, "Token", "3LED", "OmittedValues", 4.4);
        end
    end

    methods (Test)
        function exportMeasurement(testCase)
            mObj = testCase.OBJECT.RED.Measurements(end);
            mObj.exportToJson(fullfile(testCase.OUTPUT_PATH, "measurement.json"));

            jObj = pattersonlab.core.color.SpectralMeasurement.initFromJSON(...
                fullfile(testCase.OUTPUT_PATH, "measurement.json"));

            testCase.verifyEqual(mObj.wavelengths, jObj.wavelengths, "Wavelengths do not match!");
            testCase.verifyEqual(mObj.spectra0, jObj.spectra0, "Spectra do not match!");
            testCase.verifyEqual(mObj.voltage, jObj.voltage, "Voltages do not match!");
            testCase.verifyEqual(mObj.beamDiameter, jObj.beamDiameter, "Beam diameters do not match!");
            testCase.verifyEqual(mObj.minWL, jObj.minWL, "Min wavelengths do not match!");
            testCase.verifyEqual(mObj.maxWL, jObj.maxWL, "Max wavelengths do not match!");
            testCase.verifyEqual(mObj.NDF, jObj.NDF, "NDFs do not match!");
            testCase.verifyEqual(mObj.backgroundSpectra, jObj.backgroundSpectra, "Background spectra do not match!");
        end

        function exportGammaRamp(testCase)
            rObj = testCase.OBJECT.RED;
            rObj.exportToJson(fullfile(testCase.OUTPUT_PATH, "gammaramp.json"));

            jObj = pattersonlab.core.color.GammaRampMeasurement.initFromJSON(...
                fullfile(testCase.OUTPUT_PATH, "gammaramp.json"));

            testCase.verifyClass(rObj.Measurements, 'pattersonlab.core.color.SpectralMeasurement');
            testCase.verifyEqual(rObj.values, jObj.values,...
                "Mesurement values do not match!",...
                "AbsTol", 1e-5);
            testCase.verifyEqual(rObj.beamDiameter, jObj.beamDiameter,...
                "Beam diameters do not match!");
            testCase.verifyEqual(rObj.importProps.MinWL, jObj.importProps.MinWL, ...
                "Min wavelengths do not match!");
            testCase.verifyEqual(rObj.importProps.MaxWL, jObj.importProps.MaxWL, ...
                "Max wavelengths do not match!");
            testCase.verifyEqual(rObj.NDF, jObj.NDF, "NDFs do not match!");
            testCase.verifyEqual(rObj.omittedValues, jObj.omittedValues,...
                "Omitted values do not match!");

            % Check measurement order
            testCase.verifyEqual(rObj.Measurements(1).voltage, jObj.Measurements(1).voltage, ...
                "Measurement indexing does not match!");

            % Check derived properties
            testCase.verifyEqual(rObj.numMeasurements, jObj.numMeasurements,...
                "Number of measurements do not match!");
            testCase.verifyEqual(rObj.maxPower, jObj.maxPower, "Max power does not match!");
            testCase.verifyEqual(rObj.peakWavelength, jObj.peakWavelength, ...
                "Peak wavelengths do not match!");
            testCase.verifyEqual(rObj.spectra, jObj.spectra, "Calculated spectra do not match!");
            testCase.verifyEqual(rObj.wavelengths, jObj.wavelengths, "Wavelengths do not match!");
            testCase.verifyEqual(rObj.powers, jObj.powers, "Calculated powers do not match!");
        end

        function exportLedCalibration(testCase)
            cObj = testCase.OBJECT;
            cObj.exportToJSON(testCase.OUTPUT_PATH);
            jObj = pattersonlab.core.color.LedCalibration.initFromJSON(...
                fullfile(testCase.OUTPUT_PATH, "LedCalibration_20250402_0ndf.json"));

            % Check subclasses
            testCase.verifyClass(jObj.RED, 'pattersonlab.core.color.GammaRampMeasurement');
        end
    end
end