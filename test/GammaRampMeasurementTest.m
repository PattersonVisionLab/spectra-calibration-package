classdef GammaRampMeasurementTest < matlab.unittest.TestCase

    properties
        OBJECT
        DATA_PATH
    end

    methods (TestClassSetup)
        function methodSetup(testCase)

            import pattersonlab.core.color.*;

            testCase.DATA_PATH = fullfile(util.getMainFolder(),...
                "test", "data", "20250402_LEDSpectra");
            testCase.OBJECT = GammaRampMeasurement('530nm',...
                testCase.DATA_PATH, 5:-0.2:0, 3.9, "20250402", "NDF", 0,...
                'MinWL', 525, 'MaxWL', 560, "Token", "3LED", "OmittedValues", 4.4);
        end
    end

    methods (Test)
        function noMinMaxWL(testCase)
            import pattersonlab.core.color.*;
            out = GammaRampMeasurement('530nm',...
                testCase.DATA_PATH, 5:-0.2:0, 3.9, "20250402", "NDF", 0,...
                "Token", "3LED", "OmittedValues", 4.4);
        end

        function setNDF(testCase)
            import pattersonlab.core.color.*;
            obj = GammaRampMeasurement('530nm',...
                testCase.DATA_PATH, 5:-0.2:0, 3.9, "20250402", "NDF", 0,...
                 'MinWL', 525, 'MaxWL', 560, "Token", "3LED", "OmittedValues", 4.4);
            obj.setNDF(1);
            testCase.verifyEqual(obj.NDF, 1, "NDF should be set to 1!");
            testCase.verifyGreaterThan(testCase.OBJECT.maxPower, obj.maxPower,...
                "Mean power with NDF should be less than without NDF!");
        end

        function fitType(testCase)
            testCase.verifyEqual(testCase.OBJECT.lutFitType, "none");
            T1 = testCase.OBJECT.getLUT(5:-0.1:0);
            T2 = testCase.OBJECT.getLUT(5:-0.1:0, 'poly8');
            testCase.verifyNotEqual(T1.Power, T2.Power);
            testCase.verifyEqual(testCase.OBJECT.lutFitType, "poly8");
        end
    end
end