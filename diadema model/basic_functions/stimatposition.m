function [newstim] = stimatposition(stim0,stim,animalposition)
% [newstim] = stimatposition(stim0,stim,animalposition)
% This function returns the stimulus w.r.t. the animal given the position
% of animal in the arena. See section A.5 and Fig. S3 for detail.
% 
% Inputs:
%   stim0: stimulus in the arena's coordinate system (X0 in paper)
%   stim: structure, stimulus w.r.t. the arena (X in paper)
%       phi: the logitudinal location of stimuli (dots on the wall of
%       arena) w.r.t. the arena
%       intensitiy: the intensity of stimulus at each location of the arena
%       (phi)
%   animalposition: (theta, r), the polar coordinates of the final position of
%       animal (scale the radius of the arena to be 1).
% Output:
%   newstim: structure, stimulus w.r.t. the animal
%       phi: the logitudinal location of stimuli (dots on the wall of
%       arena) w.r.t. the animal
%       intensitiy: the intensity of stimulus at each phi
% All angles are with the unit of degree.
% 
% Tianshu Li
% June 15th, 2021

theta = animalposition(1); % angle of the animal's position
r = animalposition(2); % distance to the center of the arena

if r > 1
    error('Please scale the radius of the arena to 1. (animalposition(2) should not be larger than 1.)');
end

newstim = stim; % copy the original stimulus, only stim.phi changes when animal moves

original_coordinate_complex = cosd(stim0.phi)+1i*sind(stim0.phi);
animalposition_complex = r*(cosd(theta)+1i*sind(theta));
new_coordinate_complex = original_coordinate_complex-animalposition_complex;
newphi = mod(rad2deg(angle(new_coordinate_complex))+stim.center,360); % stim.center = -psi
[newphi,idx] = sort(newphi);
newphi_extended = [newphi-360,newphi,newphi+360]; % extend to negative angles and angles larger than 360 for fitting
intensity_extended = [stim0.intensity(idx), stim0.intensity(idx), stim0.intensity(idx)];
newstim.intensity = interp1(newphi_extended,intensity_extended,stim.phi);
end
