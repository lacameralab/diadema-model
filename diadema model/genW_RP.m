function [net] = genW_RP(net)
% [net] = genW_RP(net)
% This function generates the weight matrices between PRC and RN in each
% ambulacrum.
% Inputs:
%   net: structure, containing network parameters
% 
% Tianshu Li
% May 26, 2021

for kAmb = 1:net.nAmbulacrum
    if net.RN(kAmb).nCell~=net.PRC(kAmb).nCell
        error('Lateral inhibition requires same number of PRC and RN in each ambulacrum.\n');
    end
    % connection from PRC to RN
    net.param.J_RP(kAmb) = -net.param.J0_RP(kAmb); % inhibition
    net.RN(kAmb).W_RP = net.param.J_RP(kAmb)*eye(net.RN(kAmb).nCell);
    % connection between RN
    a = 0.25; % a_RR, strength of the lateral connection weights in the RNs
    net.RN(kAmb).W_RR = zeros(net.RN(kAmb).nCell,net.RN(kAmb).nCell);
    net.RN(kAmb).W_RR(1,2) = -2*a*net.param.J_RP(kAmb);
    net.RN(kAmb).W_RR(end,end-1) = -2*a*net.param.J_RP(kAmb);
    for i = 2:net.RN(kAmb).nCell-1
        net.RN(kAmb).W_RR(i,i+1) = -a*net.param.J_RP(kAmb);
        net.RN(kAmb).W_RR(i,i-1) = -a*net.param.J_RP(kAmb);
    end
end

% preferred direction of RN (similar to ONR, the direction of stimulus where RN is inhibited most. See sec. A3 for detail.)
for kAmb = 1:net.nAmbulacrum
    net.RN(kAmb).phi_pref = zeros(net.RN(kAmb).nCell,1); % preferred direction (longitudinal)
    for i = 1:net.RN(kAmb).nCell
        z = cosd(net.PRC(kAmb).phi_dms)+1i*sind(net.PRC(kAmb).phi_dms); % preferred direction as a unit vector of each PRC on the ambulacrum
        net.RN(kAmb).phi_pref(i) = mod(rad2deg(angle(sum(abs(net.RN(kAmb).W_RP(i,:))'.*z))),360);
    end
end
end