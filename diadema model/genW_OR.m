function [net] = genW_OR(net)
% [net] = genW_OR(net)
% This function generates the weight matrices between RNs in all ambulacra
% and ONR cells.
% Inputs:
%   net: structure, containing network parameters
% 
% Tianshu Li
% May 26, 2021

nRN = net.RNind(end)-1; % total number of RNs
net.ONR.W_OR = zeros(net.ONR.nCell,nRN); % synaptic weight matrix from RN to ONR

if nRN~=net.ONR.nCell
    error('Lateral inhibition requires same number of RN (all ambulacra) and ONR.\n');
end
% connection from RN to ONR
net.param.J0_OR = 1/sqrt(net.ONR.nCell); % balanced net, scaling factor of synaptic weight from RN to ONR
net.param.J_OR = -net.param.J0_OR;
net.ONR.W_OR = net.param.J_OR*eye(net.ONR.nCell);
% connection between ONR
a = 0.25; % a_OO, strength of the lateral connection weights in the ONR
net.ONR.W_OO = zeros(net.ONR.nCell);
net.ONR.W_OO(1,2) = -a*net.param.J_OR;
net.ONR.W_OO(1,end) = -a*net.param.J_OR;
net.ONR.W_OO(end,end-1) = -a*net.param.J_OR;
net.ONR.W_OO(end,1) = -a*net.param.J_OR;
for i = 2:net.ONR.nCell-1
    net.ONR.W_OO(i,i+1) = -a*net.param.J_OR;
    net.ONR.W_OO(i,i-1) = -a*net.param.J_OR;
end

% direction of maximum sensitivity of ONR (the direction of stimulus where ONR is inhibited most. See sec. A3 for detail.)
z = [];  % preferred direction vector of each RN on all ambulacra
for kAmb = 1:net.nAmbulacrum
    z = [z;cosd(net.RN(kAmb).phi_pref)+1i*sind(net.RN(kAmb).phi_pref)];
end
net.ONR.phi_pref = zeros(net.ONR.nCell,1); % direction of maximum sensitivity
for i = 1:net.ONR.nCell
    net.ONR.phi_pref(i) = mod(rad2deg(angle(sum(abs(net.ONR.W_OR(i,:))'.*z))),360); % the direction of input where ONR cell responds most (is inhibited most).
end

end