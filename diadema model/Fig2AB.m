% This script generates panel A and B of Figure 2 in the paper.
%
% Tianshu Li
% Jan. 4, 2023

clear;

addpath('./utils/');

%% figure setting
col = lines;
setcol = 0.2*ones(7,4);
setcol(:,1:3) = col(1:7,:);
fntsz = 18; % fontsize

%% Network parameters
% Structure parameters
net.nAmbulacrum = 5; % number of ambulacra
net.ONR.nCell = 500; % number of ONR cells on oral nerve ring
net.ONR.loc = linspace(0,360,net.ONR.nCell+1); % location of each ONR cell
net.ONR.loc = net.ONR.loc(1:end-1)'; % 360 is the same as 0

for kAmb = 1:net.nAmbulacrum
    % number of cells on each ambulacrum
    net.PRC(kAmb).nCell = 100; % number of PRC on each ambulacrum
    net.RN(kAmb).nCell = 100; % number of RN on each ambulacrum, bottleneck in structure (fewer neurons in higher level)
end
RNind = [1,1+cumsum([net.RN(:).nCell])]; % first index of RN in each ambulacrum
nRN = RNind(end)-1; % total number of RNs
net.RNind = RNind;

% Physiology parameters
% Distribution of the PRCs (PRCdistrtype), sigma defines the width of the distribution, which is an uniform distribution
net.param.delta = 15; % degree, Half-width of the distribution of PRCs locations on each ambulacrum
net.param.Delta_rho = 60-2*net.param.delta;%16; % degree, Acceptance angle of PRCs (width at half max of angular sensitivity function)
% Angular sensitivity curves, Delta_rho defines the width of the curve
net.param.aPRCsensitivitycurve = 2*cosd(net.param.Delta_rho/2)-1; % a, parameter controls the width of the cosine angular sensitivity curves. (-inf<a<1) When a increases, the width of the curve decreases. When a < -1, response is always positive (PRC responds to signal from all directions).

% Maximal activity of PRC, RN, and ONR
net.param.rPRCmax = 1; % r^{PRC}_{max}, maximal activity of PRCs
net.param.rRNmax = 1; % r^{RN}_{max}, maximal activity of RNs
net.param.rONRmax = 1; % r^{ONR}_{max}, maximal activity of ONRs

for kAmb = 1:net.nAmbulacrum
    % PRC, the direction of maximum sensitivity
    net.PRC(kAmb).phi_dms = linspace(-net.param.delta,net.param.delta,net.PRC(kAmb).nCell) + (kAmb-1)*360/net.nAmbulacrum; % set the location of PRCs uniform to represent the distribution in plots
    net.PRC(kAmb).phi_dms(net.PRC(kAmb).phi_dms<0) = net.PRC(kAmb).phi_dms(net.PRC(kAmb).phi_dms<0)+360; % change negative values to positive (e.g. -1 to 359)
    net.PRC(kAmb).phi_dms = net.PRC(kAmb).phi_dms'; % row vector to column vector
    % maximal activity of PRC and RN
    net.PRC(kAmb).rmax = net.param.rPRCmax*ones(net.PRC(kAmb).nCell,1); %  % maximal activity of each PRC
    net.RN(kAmb).rmax = net.param.rRNmax*ones(net.RN(kAmb).nCell,1); % maximal activity of each RN cell
    % acceptance angle of each PRC
    net.PRC(kAmb).Delta_rho = net.param.Delta_rho*ones(1,net.PRC(kAmb).nCell);
    net.PRC(kAmb).aPRCsensitivitycurve = 2*cosd(net.PRC(kAmb).Delta_rho/2)-1; % a, parameter controls the width of the tuning curve. (-inf<a<1) When a increases, the width of tuning curve decreases. When a < -1, response is always positive (PRC responds to signal from all directions).
end
net.ONR.rmax = net.param.rONRmax*ones(net.ONR.nCell,1); % maximal activity of each ONR cell

%% stimulus
acc = 0.1; % degress, bin size for stimulus
stim.phi = 0:acc:360-acc;
stim.intensity = ones(size(stim.phi));
stim.center = 0; % deg, center of stimulus

%% Fig. 2A: Angular sensitivity curves of PRC
figure(1); clf; hold on;
for kAmb = 1:net.nAmbulacrum
    for i = 1:net.PRC(kAmb).nCell
       y = max(sensitivitycurve(stim.phi,net.PRC(kAmb).phi_dms(i),net.PRC(kAmb).aPRCsensitivitycurve(i)),0);
       
       plot(stim.phi,y,'color',setcol(kAmb,:));
       if i == round(net.PRC(kAmb).nCell/2)
           plot(stim.phi,y,'color',setcol(kAmb,1:3),'linewidth',1.5);
       end
       figset(gca,'','',fntsz);
       xlim([0,360]);
       box on;
    end
end
figset(gca,'\phi','f_i^k',18);
title('f_i^k(\phi)');

%% Fig. 2B: sigmoidal function S used to model the output of RN and ONR neurons
beta = 4.5; % Steepness of RN’s sigmoidal response function
xc = -0.45; % Location parameter of ONR's sigmoidal response function
figure(2); clf; hold on;
x = -1:0.001:1;
y = (1+tanh(beta*(x-xc)))/2;
plot(x,y,'k','linewidth',2);
figset(gca,'','',18);
figset(gca,'x','S',18);
title('S(x)');

