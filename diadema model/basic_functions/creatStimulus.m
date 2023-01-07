function [stim] = creatStimulus(stim_pattern,stimwidth,stimcenter,acc)
% [stim] = creatStimulus(stim_pattern,stimwidth,stimcenter,acc)
% This function is used to creat stimulus for sea urchin vision simulation.
% Inputs:
%   stim_pattern: string, stimulus type. Details see function
%              stimulus_create
%   stimwidth: stimulus width (phi_stim, stim_pattern wavelengths in degree)
%   stimcenter: the longitudinal center of stim_pattern, unit: degree, default:
%              0 degree.
%   acc: longitudinal bin size. unit: degree, default: 0.01 degree

unitreflectance = 1; % scale stim.intensity, reflactance (r) = relative reflactance (rr) * unitreflectance

if nargin < 4
    stimcenter = 0;
end
if nargin < 5
    acc = 0.01; % longitudinal bin size, degree
end

va          = -180:acc:180-acc;     % visual angles over which the stim_pattern is observed
ratio       = 1;                % amplitude ratio of inner and outer stripes

if strcmp(stim_pattern,'bar') || strcmp(stim_pattern,'square')
    % function stimulus_create creates bar and square stimuli with bar width stimwidth/2
    intensity = stimulus_create(stim_pattern, stimwidth*2, va, ratio) + 1; % create the stimulus (the unscaled intensity for each visual angle in va)
else
    intensity = stimulus_create(stim_pattern, stimwidth, va, ratio) + 1; % create the stimulus (the unscaled intensity for each visual angle in va)
end

% Intensity is proportional to ink value (Kirwan 2018). The
% peak of intensity is 0% ink value. And minimun intensity is 100% ink
% value.
ink = (intensity-mean(intensity))/(max(intensity)-min(intensity));
ink = 1+min(ink)-ink;
rr = 1-ink; % relative reflactance, black is 0
r = rr*unitreflectance; % reflactance
r = r*(1-0.176)+0.176; % the one used for Kirwan 2018 paper

stim.phi = mod((stimcenter-180):acc:(stimcenter+180-acc),360);
[stim.phi,idx] = sort(stim.phi);
stim.type = stim_pattern;
stim.intensity = r(idx);
stim.width = stimwidth; % deg, width of stimulus
stim.center = stimcenter; % deg, center of stimulus
stim.acc = acc;
end