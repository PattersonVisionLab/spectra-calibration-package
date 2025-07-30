
% Add spectra-calibration-tools to MATLAB's search path
addpath(genpath('..\spectra-calibration-tools'));  % Replace .. with path
% Get it at https://github.com/PattersonVisionLab/spectra-calibration-tools
% or off the lab Box folder.

import pattersonlab.core.*;

measurementFolder = fullfile(color.util.getMainFolder(), 'test', 'data', '20250409_LEDSpectra');
outputFolder = fullfile(color.util.getMainFolder(), 'output');

plotFlag = true;

calDate = "20250409";
voltages = 5:-0.2:0;  % V
beamDiameter = 3.9;   % mm
ndf = 1;              % O.D. of the neutral density filter
% Set NDF to zero if using no filter or using measured values for the NDF
% transmission rather than the published values from the manufacturer.

[obj, redLED, greenLED, blueLED] = color.io.loadMaxwellianViewSpectra(...
    measurementFolder, voltages, calDate, beamDiameter, ndf, ...
    "Token", "3LED", "OmittedValues", 4.4);

if plotFlag
    redLED.plotSpectra("ShowCutoffs", true);
    greenLED.plotSpectra("ShowCutoffs", true);
    blueLED.plotSpectra("ShowCutoffs", true);
end

% Raw lookup tables
if plotFlag
    redLED.plotLUT(voltages);
    greenLED.plotLUT(voltages);
    blueLED.plotLUT(voltages);
end

% Set the fits for each lookup table
redLED.setFitType('linear');
greenLED.setFitType('poly8');
blueLED.setFitType('linear');

% Plot again to see fitted LUTs
if plotFlag
    redLED.plotLUT(voltages);
    greenLED.plotLUT(voltages);
    blueLED.plotLUT(voltages);
end

% Write the lookup tables
redLED.writeLUT(outputFolder, voltages, 'linear');
greenLED.writeLUT(outputFolder, voltages, 'poly8');
blueLED.writeLUT(outputFolder, voltages, 'linear');

% Useful mainly if you want to have the spectra handy for other programs
redLED.writeSpectra(outputFolder);
greenLED.writeSpectra(outputFolder);
blueLED.writeSpectra(outputFolder);

% Save a JSON file containing the LED calibration data
obj.exportToJSON();

% The JSON file is saved to the output folder. Include a file path as below
% to save to a different location:
% obj.exportToJSON('C:\Users\YourName\Documents');

% You can reload the LedCalibration object from the JSON file using:
jsonObj = LedCalibration.initFromJSON(fullfile(outputFolder, "LedCalibration_20250409_10ndf.json"));


%% Process a stimulus
stim = [0.5+zeros(size(0.002:0.002:1)), sin(4*2*pi*(0:0.002:1))*0.5+0.5];
figure(); plot(0:0.002:2, stim); title('Normalized stimulus');

% This converts it into RGB modulations around a white point
ledValues = obj.calcStimulus(stim);  % uW
% Every time you convert, a plot of the RGB values is generated to sanity
% check the result before writing it to a stimulus file for an experiment.

fName = 'mystim.txt';
color.io.makeLEDStimulusFile(fName, ledValues, obj);
% The stimulus file is saved in the output folder. Include a file path ahead
% of the fName to save elsewhere.