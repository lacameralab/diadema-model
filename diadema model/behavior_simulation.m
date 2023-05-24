% This script simulate the behavior expeiment. It generates Fig. 2C and
% Fig. 4B.
%
% P.S.: the units of all angles here are "degree" ([0, 360)).
%
% Tianshu Li
% May. 25, 2021

clear;

addpath('./utils/');

%% set stimulus
% Select the type (stim_pattern) of the stimulus (please refer to the mail paper for the description of the stimuli)
% stim_pattern = 'bar'; stim.width = 40;  % deg, width of stimulus
stim_pattern = 'DoG'; stim.width = 69;  % 29; % deg, width of stimulus
% stim_pattern = 'square'; stim.width = 69;  % deg, width of stimulus
% stim_pattern = 'Hermitian'; stim.width = 69;  % deg, width of stimulus
% stim_pattern = '2square'; stim.width = 69;  % deg, width of stimulus
% stim_pattern = 'Morlet'; stim.width = 69;  % deg, width of stimulus
% stim_pattern = 'control';

% create stimulus
acc = 0.1; % degress, bin size for stimulus
stim.phi = 0:acc:360-acc;
stim.intensity = ones(size(stim.phi));
stim.center = 0; % deg, center of stimulus

if strcmp(stim_pattern,'control')
    stim = createStimulus('bar',10,stim.center,acc);
    stim.intensity = 0.77*ones(size(stim.phi)); % control stimulus
else
    stim = createStimulus(stim_pattern,stim.width,stim.center,acc);
end

%% figure setting
col = lines;
setcol = 0.2*ones(7,4);
setcol(:,1:3) = col(1:7,:);
fntsz = 18; % fontsize

%% ================= set network parameters  ==================
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
    net.PRC(kAmb).phi_dms = net.param.delta*(2*rand(net.PRC(kAmb).nCell,1)-1) + (kAmb-1)*360/net.nAmbulacrum;
    net.PRC(kAmb).phi_dms = sort(net.PRC(kAmb).phi_dms);
    net.PRC(kAmb).phi_dms(net.PRC(kAmb).phi_dms<0) = net.PRC(kAmb).phi_dms(net.PRC(kAmb).phi_dms<0)+360; % change negative values to positive (e.g. -1 to 359)
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

%% =================          simulate behavior         ===================
% One session includes Ntrials trials. The animal is placed at random
% orientation in the center of the arena in each trial, which was the case
% in experiment.
% See section A.2 and Fig. S1 for detail.
Ntrials = 100; % number of animals
stimcenters = 360*rand([1,Ntrials]); % psi, the center of stimulus with respect to the animal
vpop_amplitude = zeros(Ntrials,1); % |v_pop|, amplitude of the population vector
vpop_direction = vpop_amplitude; % <v_pop>, direction of the population vector
theta_p = 5; % if the length of vpop (population vector) is less than this threshold, the animal moves randomly

for ktrial = 1:Ntrials
    stim.center = stimcenters(ktrial); % deg, center of stimulus
    if strcmp(stim_pattern,'control')
        stim.intensity = 0.72*ones(size(stim.phi)); % control stimulus
    else
        stim = createStimulus(stim_pattern,stim.width,stim.center,acc);
    end

    net = PRCoutput(net,stim);
    [net] = RNoutput(net);
    [net.ONR] = ONRoutput(net);

    vpop = sum((net.ONR.out).*(cosd(net.ONR.phi_pref) + sind(net.ONR.phi_pref)*1i))/sqrt(net.ONR.nCell); % estimated direction of light by the animal, vpop is written as a complex number to represent a vector in Cartesian space.
    vpop_amplitude(ktrial) = abs(vpop);
    vpop_direction(ktrial) = mod(-stimcenters(ktrial)+rad2deg(angle(vpop)),360); % movement direction
end

%% population vectors (Fig. 2C)
figure(1); clf;
stim = createStimulus(stim_pattern,stim.width,0,acc);
% threshold
h2 = polarplot(linspace(0,2*pi,500),theta_p*ones(1,500),'linewidth',2);
hold on;

% population vector
h1 = polarscatter(deg2rad(vpop_direction),vpop_amplitude,10,'linewidth',2);
rmax = rlim;
rmax = rmax(2);

% stimulus
for i = 1:10
    [strike,dip]= meshgrid(stim.phi*pi/180,9+(i-1)*0.1);
    polarscatter(strike,dip,50*ones(size(strike)),repmat(stim.intensity',1,3),'filled','MarkerFaceAlpha',0.8);
end
polarplot(stim.phi*pi/180,(9-0.2)*ones(size(stim.phi)),'k','linewidth',2);
polarplot(stim.phi*pi/180,(10)*ones(size(stim.phi)),'k','linewidth',2);
colormap bone;

rlim([0,10])
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'counterclockwise';
pax.FontSize = fntsz;

title(sprintf('%i trials, %i deg %s, Delta rho = %i deg, delta = %i deg',Ntrials,stim.width,stim_pattern,net.PRC(1).Delta_rho(1),net.param.delta));

legend([h1,h2],'v_{pop}','\theta_p');
legend 'boxoff';

%% plot final position (animal heading) in all trials (Fig. 4B)
% See section 5.3.5 for details.
finalP = zeros(1,Ntrials); % final position of each animal
for i = 1:Ntrials
    if vpop_amplitude(i)<=theta_p
        finalP(i) = 360*rand(); % below threshold: sampled from a circularly uniform distribution
    else
        finalP(i) = vpop_direction(i)+randn()/(vpop_amplitude(i)-theta_p); % above threshold: sample from a Gaussian distribution narrowly centered around the population vector with standard deviation 1/(|vpop|−theta_p)
    end
end
finalP = mod(finalP, 360);
finalP = deg2rad(finalP);

figure(2); clf;
% stimulus
stim = createStimulus(stim_pattern,stim.width,0,acc);
if strcmp(stim_pattern,'control')
    stim.intensity  = 0.72*ones(size(stim.intensity)); % control
end
for i = 1:10
    [strike,dip]= meshgrid(stim.phi*pi/180,1+(i-1)*0.01);
    polarscatter(strike,dip,50*ones(size(strike)),repmat(stim.intensity',1,3),'filled','MarkerFaceAlpha',0.8);
    hold on;
end
polarplot(stim.phi*pi/180,(1-0.02)*ones(size(stim.phi)),'k','linewidth',2);
polarplot(stim.phi*pi/180,(1.12)*ones(size(stim.phi)),'k','linewidth',2);
colormap bone;

pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'counterclockwise';
pax.FontSize = 22;
rticklabels('');
rlim([0,1.12]);

hold on;
% plot final position
edgesP = linspace(0,2*pi,60); % use this vector as a seive to seperate the edge of area into small sections, stack the animals in each section to make the plot in experiment paper
finalP_sort = sort(finalP);
stackspace = 0.08; % the r difference in plot between animals in different row (stacking)
Ntrials_in_section = zeros(1,length(edgesP)-1); % number of animals in each section
while finalP_sort
    for i = 1:length(edgesP)-1
        if finalP_sort(1) >= edgesP(i) && finalP_sort(1) < edgesP(i+1)
            h1 = polarscatter(finalP_sort(1),9.3/10-stackspace*Ntrials_in_section(i),200,'b','filled','MarkerFaceAlpha',.5);
            finalP_sort(1) = [];
            Ntrials_in_section(i) = Ntrials_in_section(i)+1;
        end
        if isempty(finalP_sort)
            break;
        end
    end
end

% plot the arrow at the center, representing the statistical
% direction of the population
finalPv = [finalP', ones(length(finalP),1)]; % finalP as a vector, first column is the direction from the center to the anmial, second column is 1
[finalPvcartr, finalPvcarti] = pol2cart(finalPv(:,1),finalPv(:,2)); % change polar to cartisian coordinate
finalPvcart = finalPvcartr + finalPvcarti*1i; % change to complex number to do the calculation
popfinalP = sum(finalPvcart);
h2 = polarplot([0,angle(popfinalP)],[0,norm(popfinalP)/Ntrials],'r','linewidth',5);

legend([h1,h2],'animal heading','mean heading arrow');
legend 'boxoff';
title(sprintf('%i trials, %i deg %s, Delta rho = %i deg, delta = %i deg',Ntrials,stim.width,stim_pattern,net.PRC(1).Delta_rho(1),net.param.delta));
