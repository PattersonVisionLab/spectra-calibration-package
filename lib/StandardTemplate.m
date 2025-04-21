function spectra = StandardTemplate(lamda_max, wavelengths)
%{
    Standard template formula for computing photoreceptor sensitivity of A1 (inc. 11-cis-retinal) pigments.
    The equations below were derived originally from Lamb (1995) then later modified by
    Govardovskii etal 2000 Vis Neruosci. Govardovskii made the modifications based on fits to opsin 
    absorbtion of many animals (not including mammals) using microspectrophotometry. The original 
    equation was derived from bovine opsin.

    Input:
        lamda_max   = lamda mawavelengths of the opsin (nm, scalar)
        wavelengths = [optional] list of wavelengths to compute the spectra over (nm, column vector) 

    list of kown lamda max:
        mouse mcone = 508nm


%}

    % input handeling
        if nargin == 0
            error('Not enough arguments')
        end
        if nargin == 1 && ~exist('wavelengths', 'var')
            wavelengths = (350 : 1: 700)';
        end
    
    % Mansfield?MacNichol (MM) transform
        wavelengths_mm = lamda_max./wavelengths; 

    % construct alpha band of curve (long wavelength side)   
        % constants
            A = 67.7;
            % a = 0.88;
            a = 0.8795 + 0.0459 * exp(-(lamda_max - 3000)^2 / 11940); % modification by govardoskii etal 2000
            B = 28;
            b = 0.922;
            C = -14.9;
            c = 1.104;
            D = 0.674;

        % standard pigment absorbance template from lamb (1995)
        y = 1./( exp(A.*(a-wavelengths_mm)) + exp(B.*(b-wavelengths_mm)) + exp(C.*(c-wavelengths_mm)) + D );

    % construct beta band of curve (short wavelength secondary peak)    
        A_beta = 0.26;
        lamda_max_beta = 189 + 0.315 * lamda_max;
        beta = -40.5 + 0.195 * lamda_max;

        y_beta = A_beta .* exp(-((wavelengths - lamda_max_beta)./ beta).^2 ); 

    % complete template
        spectra = y + y_beta; % complete template of the photoreceptor sensitivity (absorbance)

end