% This script generates Figure S2 in the paper.
%
% Tianshu Li
% Jan. 4, 2023
 
clear;

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

% Connectivity parameters
% scaling factor
net.param.J0_RP = zeros(1,net.nAmbulacrum); % scaling factor of synaptic weight from PRC to RN
for kAmb = 1:net.nAmbulacrum
    net.param.J0_RP(kAmb) = 1/sqrt(net.PRC(kAmb).nCell);
end
net.param.J0_OR = 1/sqrt(net.nAmbulacrum*net.RN(1).nCell); % scaling factor of synaptic weight from RN to ONR

% generate W_RP
net = genW_RP(net);
% generate W_OR
net = genW_OR(net);

%% Response function parameters
% RN output
beta_RN = 3; % Steepness of RN’s sigmoidal response function
xc_RN = -0.6; % % Location parameter of RN’s sigmoidal response function
net.param.beta_RN = beta_RN;
net.param.xc_RN = xc_RN;
% ONR output
beta_ONR = 4.5; % Steepness of RN’s sigmoidal response function
xc_ONR = -0.45; % Location parameter of ONR's sigmoidal response function
net.param.beta_ONR = beta_ONR;
net.param.xc_ONR = xc_ONR;

%% calculate preferred directions of eONR cells by simulation
acc = 0.1; % degress, bin size for stimulus

% Preferred direction is defined as the direction of a narrow stimulus
% causing the maximal increase in activity in the cell. Since eONR are
% inhibited by light, the preferred direction is the direction of a narrow
% light stimulus causing the maximal decrease in activity in the eONR cell.
% Use a very narrow bright strip with intensity 1 at strip and 0 elsewhere
% for the narrow stimulus.
% To find the preferred direction of an ONR cell:
% 1. Simulate the output of this ONR cell given this kind of stimulus
% centering at all directions.
% 2. Find the stimulus that inhibits this cell most.

% create the record for each cell
ONR = net.ONR;
for i = 1:net.ONR.nCell
    ONR(i).stimcenter = []; % record the center of the stimulus (direction of the dot light)
    ONR(i).out = []; % record the corresponding output of cell i
end
% simulation
for stimcenter = 1:360
    stim.width = 2; % deg, width of the white strip
    stim.center = stimcenter; % deg, center of stimulus
    stim = creatStimulus('bar',360-stim.width,mod(180+stim.center,360),acc); % use 360-stim.width bar to generate white strip
    stim.center = stimcenter; % deg, center of stimulus
    stim.intensity = stim.intensity<0.5; % set the light intensity of the black area to 0

    % simulate the output of all ONR cells
    net = PRCoutput(net,stim);
    net = RNoutput(net);
    net.ONR = ONRoutput(net);

    % add record to each cell
    for kONR = 1:net.ONR.nCell
        ONR(kONR).stimcenter = [ONR(kONR).stimcenter,stimcenter]; % record the center of the stimulus (direction of the dot light)
        ONR(kONR).out = [ONR(kONR).out,net.ONR.out(kONR)]; % record the corresponding output of cell kONR
    end
end

net.ONR.phi_pref_sim = zeros(net.ONR.nCell,1); % preferred direction calculated by simulation
for kONR = 1:net.ONR.nCell
    [y,idx] = min(ONR(kONR).out); % find the direction where the output of cell kONE is the lowest
    x = ONR(kONR).stimcenter(idx);
    net.ONR.phi_pref_sim(kONR) = x;
end

%% simulated vs analytical preferred direction of ONRs
figure(3); clf;
hold on;
% simulated preferred direction
plot(net.ONR.loc, net.ONR.phi_pref_sim,'x','markersize',10);
% analytical preferred direction (calculated in function genW_OR)
plot(net.ONR.loc, net.ONR.phi_pref,'.','markersize',10);
xticks(0:72:360);
yticks(0:72:360);
xlim([0,360]);
ylim([0,360]);
figset(gca,'Location of eONR neurons on the oral nerve ring','preferred direction of eONR neurons',fntsz);
legend('simulated','analytical','fontsize',18);
legend 'boxoff';