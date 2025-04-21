# Color calibration package

Calibration code for 3-primary stimulators. Designed for the Maxwellian View in the 1P AOSLO.

### Layout
- __data__: transmission spectra for commonly-used ThorLabs ND filters ([source](https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=266&pn=NE03A#811))
- __output__: spectra and gamma files created (not under version control)
- __lib__: 3rd party code
- __scripts__: sample calibration routines
- __src__: main code
    - __+pattersonlab/+core/+color__: main calibration classes
            - __+io__: import/export functions
            - __+util__: support functions
            - `SpectralMeasurement`: a single measurement from the spectrometer
            - `GammaRampMeasurement`: a group of measurements with varying input values
            - `LedCalibration`: calibrations for 3 LEDs
- __test__: unit test suite

The current code coverage of the test suite is 83.82% of functions and 71.09% of statements.

### Features
- Reads measurements from the Ocean Optics STS-VIS spectrometer (`SpectralMeasurement`)
- Generate a lookup table relating voltage to power in microwatts (`GammaRampMeasurement`).
- Identifies a set of scalar weights (one per LED) needed to achieve a white point defined as 0.33 0.33 on CIE diagram (`LedCalibration`)
- Calculates change in power for 3 LEDs from the white point needed to generate achromatic, isoluminant, and cone-isolating stimuli (`LedCalibration`)
- Writes all the measurements and calibration calculation output to a single JSON file which can be used to recreate the MATLAB objects
- Creates LED stimulus files in the correct format for the one-photon AOSLO (`makeLEDStimulusFile`).

The routines in `LedCalibration` assume macaque cone peak sensitivities and macular pigment. The cone peak sensitivities can be changed when initializing `LedCalibration`, but currently macular pigment spectral data is only present for macaque.

__See also:__ Maxwellian View calibration protocols in [PattersonLabDocs](https://github.com/PattersonVisionLab/PattersonLabDocs)


### Dependencies
Requirements: MATLAB with the Curve Fitting Toolbox and the Optimization Toolbox. Data import is designed for Ocean Optics spectrometers.

3rd party code is from  [Psychtoolbox](https://github.com/Psychtoolbox-3/Psychtoolbox-3), [SilentSubstitutionToolbox](https://github.com/spitschan/SilentSubstitutionToolbox), the Neitz lab, and my other repositories.


### Citations
- Godat T, Kohout K, Parkins K, Yang Q, McGregor JE, Merigan WH, Williams DR, Patterson SS (2024) Cone-Opponent Ganglion Cells in the Primate Fovea Tuned to Noncardinal Color Directions. _Journal of Neuroscience_, 44(18), e1738232024
- Spitschan M, Aguirre GK, Brainard DH (2015) Selective Stimulation of Penumbral Cones Reveals Perception in the Shadow of Retinal Blood Vessels. _PLoS ONE_, 10(4), e0124328.