function [net] = PRCoutput(net,stim)
% [net] = PRCoutput(net,stim)
% This function calculates the output of all PRCs given a network and
% stimulus.
% Inputs:
%   net: structure, urchin network
%   stim: structure, stimulus, with field, 'phi' (light input
%       angle),'intensity' (light intensity at each angle phi), 
% 
% Tianshu Li
% May 25, 2021

if ~isfield(stim,'acc')
    stim.acc = stim.phi(2)-stim.phi(1); % degrees/pixel
end
for kAmb = 1:net.nAmbulacrum
    net.PRC(kAmb).out = zeros(net.PRC(kAmb).nCell,1);
    for i = 1:net.PRC(kAmb).nCell
        a = net.PRC(kAmb).aPRCsensitivitycurve(i);
        R = 2*(-1+a^2+a*sqrt(1-a^2)*acos(a))/((a-1)*sqrt(1-a^2)); % for normalization, integrate (angular sensitivity curve) ONLY FOR a > -1
        y = sensitivitycurve(stim.phi,net.PRC(kAmb).phi_dms(i),a)*stim.intensity'/R;
        net.PRC(kAmb).out(i) = net.PRC(kAmb).rmax(i)*(stim.acc*pi/180)*y;
    end
end
end