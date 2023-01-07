function outint = stimulus_create(type, wl, va, ratio)
% STIMULUS_CREATE generates a 2D stimulus vector for phototaxis experiments of desired type and resolution
%   outint = stimulus_create(type, wl, va, ratio)
%
% Inputs:
%   type   - type of stimulus, e.g. 'bar'/'dog'/'square'/'log' (see below for full list)
%   wl     - "width" or "wavelength" of pattern in degrees (definition varies with stimulus type)
%   va     - visual angles over which to create the stimulus (in degrees)
%   ratio  - ratio between black and white areas; usually ratio of amplitudes, but definition depends on stimulus type
%
% Outputs:
%   outint - stimulus intensity over all angles in va (arbitraty units, normalise later)
%
% Pattern type options:
%                'bar'          % a single black bar of width wl/2 on a white background
%                'singlegauss'  % single dark gaussian (with half-width wl/2) on white background
%                'square'       % square wave pattern: black stripe of width wl/2 flanked by two white stripes of half the width
%                '2square'      % square wavelet pattern (Haar wavelet): black stripe and white stripe, each of width wl/2, on grey background
%                'dog'          % difference of Gaussians, central Gaussian with half-width wl/2
%                'dog2'         % difference of Gaussians with different half-width ratio
%                'log'          % Laplacian of Gaussian, v1
%                'mylog'        % Laplacian of Gaussian, v2
%                'cos'          % piece-wise sine with wavelenth wl (as used for Onychophorans); ratio ignored, always 0.5/1/0.5
%                'sin'          % continous sine-wave pattern with wavelenth wl
%                'multiple'     % 7 progressively smaller cosines
%
% example usage: outint = stimulus_create('dog', 10, -180:0.1:180, 1)
%
% Written by Jochen Smolka and John D. Kirwan, Lund University, 2017
%
% This file is published as supplementary to the following article:
%   John D. Kirwan, Michael J. Bok, Jochen Smolka, James J. Foster, José Carlos Hernández, Dan-Eric Nilsson (2017) 
%   The sea urchin Diadema africanum uses low resolution vision to find shelter and deter enemies. Journal of Experimental Biology.

%% set defaults
if nargin<4, ratio = 1; end
if nargin<3, va = -180:0.1:180; end
if nargin<2, wl = 10; end
if nargin<1, type = 'bar'; end
    
outint = zeros(size(va)); % preallocate vector of pattern intensity

%% for all piece-wise defined patterns (e.g. square, triangle, sin), define the regions of black and white "stripes"
sel_whiteleft       = va>=-(ratio+1)*wl/2 & va< -wl/2; 
sel_black           = va>=-wl/2 & va< wl/2;
sel_whiteright      = va>= wl/2 & va< (ratio+1)*wl/2;
amps                = [1 1/ratio]; % relative amplitudes of black and white areas

%% main switch to create patterns
switch lower(type)
    %% single bar patterns
    case 'bar'          % single black bar (of width wl/2) on white background
        sel_black           = va>=-wl/4 & va< wl/4;
        outint(sel_black) = -1;
        
    case 'singlegauss'  % single dark gaussian (with half-width wl/2) on white background
        sigma  = wl/2 / (2*sqrt(2*log(2)));   % sigma of Gaussian
        outint = - exp(-va.^2/(2*sigma^2));
               
    %% balanced patterns (integral = 0)
    case 'square'       % square wave pattern: black stripe of width wl/2 flanked by two white stripes of half the width
        sel_whiteleft       = va>=-(ratio+1)*wl/4 & va< -wl/4; 
        sel_black           = va>=-wl/4 & va< wl/4;
        sel_whiteright      = va>= wl/4 & va< (ratio+1)*wl/4;
        outint(sel_black)        = -amps(1);
        outint(sel_whiteleft)    = amps(2);
        outint(sel_whiteright)   = amps(2);
        
    case '2square'      % square wavelet pattern (Haar wavelet): black stripe and white stripe, each of width wl/2, on grey background
        sel_white               = va>=-wl/2 & va< 0; 
        sel_black               = va>=0 & va< wl/2;
        outint(sel_black)       = -amps(1);
        outint(sel_white)       = amps(2);
        
    case 'dog'          % difference of Gaussians pattern
        fwhm1   = wl/2;                           % half-width of primary (black) Gaussian; also equals half the distance between white peaks
        fwhm2   = wl/2 * (ratio+1);               % half-width of secondary (white) Gaussian
        sigma1  = fwhm1 / (2*sqrt(2*log(2)));     % sigma of primary Gaussian
        sigma2  = fwhm2 / (2*sqrt(2*log(2)));     % sigma of secondary Gaussian
        amp1    = -1;
        amp2    = -amp1 / (fwhm2/fwhm1);
        outint  = amp1 * exp(-va.^2/(2*sigma1^2)) + amp2 * exp(-va.^2/(2*sigma2^2));        
    
    case 'dog2'         % difference of Gaussians pattern with different half-width ratio
        fwhm1   = wl;                           % half-width of primary (black) Gaussian; also equals half the distance between white peaks
        fwhm2   = wl * 1.2;                     % half-width of secondary (white) Gaussian
        sigma1  = fwhm1 / (2*sqrt(2*log(2)));   % sigma of primary Gaussian
        sigma2  = fwhm2 / (2*sqrt(2*log(2)));   % sigma of secondary Gaussian
        amp1    = -1*6;
        amp2    = -amp1 / (fwhm2/fwhm1);
        outint  = amp1 * exp(-va.^2/(2*sigma1^2)) + amp2 * exp(-va.^2/(2*sigma2^2));       
        
    case 'log'          % Laplacian of Gaussian (not working properly?)
        sigma   = wl/2 / (2*sqrt(2*log(2)));      % sigma
        acc     = median(diff(va));             % degrees/pixel
        outint  = fspecial('log', size(va), sigma/acc);
        
    case 'mylog'        % Laplacian of Gaussian
        sigma   = wl/2 / (2*sqrt(2*log(2)));
        outint = -(1-va.^2/sigma^2).*exp(-va.^2/(2*sigma^2));
        
    case 'cos'          % piece-wise sine with wavelenth wl (as used for Onychophorans)
        sel_whiteleft       = va>=-3/4*wl & va< -1/4*wl; 
        sel_black           = va>=-1/4*wl & va<  1/4*wl;
        sel_whiteright      = va>= 1/4*wl & va<  3/4*wl;
        amps                = [1 1/2]; % relative amplitudes of black and white areas
        outint(sel_black)        = - amps(1) * cosd(linspace(-90,90,nnz(sel_black)));
        outint(sel_whiteleft)    =   amps(2) * cosd(linspace(-90,90,nnz(sel_whiteleft)));
        outint(sel_whiteright)   =   amps(2) * cosd(linspace(-90,90,nnz(sel_whiteright)));
    
    case 'sin'          % continous sine-wave pattern with wavelenth wl
        outint = cosd(va*360/wl);
        
    case 'multiple'     % piece-wise cosine with 7 pieces (was intended to more closely simulate a continuous sine-wave)
        sel_white2left           = va>=-(3*ratio+1)*wl/2 & va<= -(2*ratio+1)*wl/2;
        sel_black2left           = va>=-(2*ratio+1)*wl/2 & va<= -(ratio+1)*wl/2;
        sel_black2right          = va>=(1*ratio+1)*wl/2 & va<= (2*ratio+1)*wl/2;
        sel_white2right          = va>=(2*ratio+1)*wl/2 & va<= (3*ratio+1)*wl/2;
        amps                     = [1 4/3/ratio 2/3/ratio 2/6/ratio];
        outint(sel_black)        = - amps(1) * cosd(linspace(-90,90,nnz(sel_black)));
        outint(sel_whiteleft)    =   amps(2) * cosd(linspace(-90,90,nnz(sel_whiteleft)));
        outint(sel_whiteright)   =   amps(2) * cosd(linspace(-90,90,nnz(sel_whiteright)));
        outint(sel_black2left)   = - amps(3) * cosd(linspace(-90,90,nnz(sel_black2left)));
        outint(sel_black2right)  = - amps(3) * cosd(linspace(-90,90,nnz(sel_black2right)));
        outint(sel_white2left)   =   amps(4) * cosd(linspace(-90,90,nnz(sel_white2left)));
        outint(sel_white2right)  =   amps(4) * cosd(linspace(-90,90,nnz(sel_white2right)));
        
    % The following is added by Tianshu Li in June 2021
    case 'hermitian'
        D = 2*wl/(2*sqrt(3)*sqrt(2*log(2)));
        Gwv = exp(-va.^2/(2*D^2));
        Gwv2 = [diff(Gwv)./diff(va),0];
        outint = (Gwv2-min(Gwv2))/(max(Gwv2)-min(Gwv2));
        
    case 'morlet'
        D = 2*wl/(2*sqrt(2*log(2)));
        phi = exp(-va.^2/(2*D^2)).*cosd(360*va/(wl));
        outint = (-phi+1)/1.8852236;
end

end %main