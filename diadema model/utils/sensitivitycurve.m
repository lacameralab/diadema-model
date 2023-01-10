function [y] = sensitivitycurve(phi, pref_phi, a)
% [y] = cosinetuning(phi, pref_phi, a)
% This function is the tuning curve of PRC (longitudinal). We use cosine
% tuning curve.
% Inputs:
%   phi: relative direction of light (signal)
%   pref_phi: preferred direction of the cell
%   a: parameter controls the width of the angular sensitivity curve.
%    (-inf<a<1) When a increases, the width of sensitivity curve decreases.
%    When a < -1, response y is always positive (PRC responds to signal from
%    all directions).
%    Default: 0.95
% 
% Reference: Salinas E, Abbott LF (1994) Vector reconstruction from firing
%   rates. J Comput Neurosci 1: 89–107. 
% 
% Tianshu Li
% Fen. 25, 2021

if nargin < 3
    a = 0.95;
end

y = (cosd(phi-pref_phi)-a)/(1-a);
y(y<0) = 0;

end