
import pattersonlab.core.color.*;

measurementFolder = fullfile(util.getMainFolder(), 'test', 'data', '20250409_LEDSpectra');
outputFolder = fullfile(util.getMainFolder(), 'output');

plotFlag = true;

calDate = "20250409";
voltages = 5:-0.2:0;  % V
beamDiameter = 3.9;   % mm
ndf = 1;              % O.D. of the neutral density filter
% Set NDF to zero if using no filter or using measured values for the NDF
% transmission rather than the published values from the manufacturer.

[obj, redLED, greenLED, blueLED] = io.loadMaxwellianViewSpectra(...
    measurementFolder, voltages, calDate, beamDiameter, ndf, ...
    "Token", "3LED", "OmittedValues", 4.4);

if plotFlag
    redLED.plotSpectra("ShowCutoffs", true);
    greenLED.plotSpectra("ShowCutoffs", true);
    blueLED.plotSpectra("ShowCutoffs", true);
end

if plotFlag
    redLED.plotLUT(voltages);
    greenLED.plotLUT(voltages);
    blueLED.plotLUT(voltages);
end

redLED.writeLUT(outputFolder, voltages, 'linear');
greenLED.writeLUT(outputFolder, voltages, 'poly8');
blueLED.writeLUT(outputFolder, voltages, 'linear');

if plotFlag
    redLED.plotLUT(voltages);
    greenLED.plotLUT(voltages);
    blueLED.plotLUT(voltages);
end

redLED.writeSpectra(outputFolder);
greenLED.writeSpectra(outputFolder);
blueLED.writeSpectra(outputFolder);

% Save a JSON file containing the LED calibration data
obj.exportToJSON();
% The JSON file is saved to the output folder. Include a file path as below
% to save to a different location:
% obj.exportToJSON('C:\Users\YourName\Documents');
% You can reload the LedCalibration object from the JSON file using:
% jsonObj = LedCalibration.initFromJSON('C:\Users\YourName\Documents\LEDCalibration.json');



% Process a stimulus
stim = [0.5+zeros(size(0.002:0.002:1)), sin(4*2*pi*(0:0.002:1))*0.5+0.5];
figure(); plot(0:0.002:2, stim); title('Normalized stimulus');

ledValues = obj.calcStimulus(stim);
fName = 'mystim.txt';
io.makeLEDStimulusFile(fName, ledValues, obj);
% The stimulus file is saved in the output folder. Include a file path ahead
% of the fName to save elsewhere.