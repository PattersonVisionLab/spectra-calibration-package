classdef LedCalibrationTest < matlab.unittest.TestCase

    properties
        OBJECT
        DATA_PATH
    end

    methods (TestClassSetup)
        function methodSetup(testCase)

            testCase.DATA_PATH = fullfile(pattersonlab.core.color.util.getMainFolder(),...
                "test", "data", "20250402_LEDSpectra");
            testCase.OBJECT = pattersonlab.core.color.io.loadMaxwellianViewSpectra(...
                testCase.DATA_PATH, 5:-0.2:0, "20250402", 3.9, 0, "Token", "3LED", "OmittedValues", 4.4);
        end
    end

    methods (Test)
        function ndfTest(testCase)
            withNDF = pattersonlab.core.color.io.loadMaxwellianViewSpectra(...
                testCase.DATA_PATH, 5:-0.2:0, "20250402", 3.9, 1, "Token", "3LED", "OmittedValues", 4.4);
            disp('Running test')

            testCase.verifyGreaterThan(testCase.OBJECT.meanPower, withNDF.meanPower,...
                "Mean power with NDF should be less than without NDF!");
        end

        function ciePlotTest(testCase)
            ax = testCase.OBJECT.plotCIE();
            delete(ax.Parent);
        end
    end
end
