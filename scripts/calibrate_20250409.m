
import pattersonlab.core.color.*;

measurementFolder = fullfile(util.getMainFolder(), 'test', 'data', '20250409_LEDSpectra');
outputFolder = fullfile(util.getMainFolder(), 'output');

plotFlag = true;


calDate = "20250409";
voltages = 5:-0.2:0;  % V
beamDiameter = 3.9;   % mm
ndf = 1;

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



% Process a stimulus
ledValues = obj.calcStimulus("Lum", stim);
io.writeLedStimulusFile();