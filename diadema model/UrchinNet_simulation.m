% This script creates a model of A SLICE OF sea urchin visual system and
% simulate one trial of the experiment.
% 5 ambulacra, PRC, RN on each ambulacrum, ONR cells on the center ring
% light --| photoreceptor cells (PRC) --| radial nerves (RN) --> oral nerve
% ring layer (ONR)
% 
% Since we do not have latitude difference in experient, we only model PRC
% and RN on a slice on horizontal axis of the animal with ONR.
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

%% visualize the stimulus
figure(1); clf;
plot(stim.phi, stim.intensity,'k','linewidth',2);
figset(gca,'\phi','stimulus intensitiy',fntsz);
xlim([0,360]);
title(sprintf('Stimulus: %g deg %s',stim.width,stim_pattern));

%% Simulation
% PRC output
[net] = PRCoutput(net,stim);
% RN output
net = RNoutput(net);
% ONR output
net.ONR = ONRoutput(net);

%% ================      output visualization       =================== 
figure(2); clf; 
set(gcf,'position',[20          73        1661         843]);

%% plot PRC responses to stimulus
subplot(2,2,1);
rmax = net.param.rPRCmax/2; % maximum r tick
% PRC responses
for kAmb = 1:net.nAmbulacrum
    h(kAmb) = polarplot(net.PRC(kAmb).phi_dms*pi/180,net.PRC(kAmb).out+net.param.rPRCmax/2,'o');
    hold on;
    rmax = max(rmax, net.param.rPRCmax/2+max(net.PRC(kAmb).out));
end

[ah] = ploturchin(gca, stim, net.param.rPRCmax/2, rmax+net.param.rPRCmax/2, fntsz);

rlim([0,rmax+net.param.rRNmax/2]);
rticks(rmax);
rticklabels(sprintf('%g',rmax-net.param.rRNmax/2));
legend(h,'r^{1,PRC}','r^{2,PRC}','r^{3,PRC}','r^{4,PRC}','r^{5,PRC}');
legend 'boxoff';
title('PRC response');
clear h;

%% plot RN responses to stimulus
subplot(2,2,2); cla; 
rmax = net.param.rRNmax/2; % maximum r tick
% RN responses
for kAmb = 1:net.nAmbulacrum
    h(kAmb) = polarscatter(net.RN(kAmb).phi_pref*pi/180,net.RN(kAmb).out+net.param.rRNmax/2,'o');
    hold on;
    rmax = max(rmax, net.param.rRNmax/2+max(net.RN(kAmb).out));
end

[ah] = ploturchin(gca, stim, net.param.rRNmax/2, rmax+net.param.rRNmax/2, fntsz);

rlim([0,rmax+net.param.rRNmax/2]);
legend(h,'r^{1,RN}','r^{2,RN}','r^{3,RN}','r^{4,RN}','r^{5,RN}');
legend 'boxoff';
title('RN response');
clear h;

%% plot ONR responses to stimulus
subplot(2,2,3); cla; 
rmax = net.param.rONRmax/2; % maximum r tick

% ONR responses
h = polarplot(net.ONR.phi_pref*pi/180,net.ONR.out+net.param.rONRmax/2,'o');
rmax = max(rmax, net.param.rONRmax/2+max(net.ONR.out));
hold on;
[ah] = ploturchin(gca, stim, net.param.rONRmax/2, rmax+net.param.rONRmax/2, fntsz);

rlim([0,rmax+net.param.rRNmax/2]);
legend(h,'r^{ONR}');
legend 'boxoff';
title('ONR response');
clear h;

%% population vector
theta_p = 5; % Threshold for population vector length, if the length of z (population vector) is less than this threshold, the animal moves randomly
pop_vect = zeros(1,2); % first column is amplitude of population vector, second column is the direction

z = sum((net.ONR.out).*(cosd(net.ONR.phi_pref) + sind(net.ONR.phi_pref)*1i))/sqrt(net.ONR.nCell); % estimated direction of light by the animal
pop_vect(1) = abs(z); % amplitude
pop_vect(2) = mod(rad2deg(angle(z)),360); % move towards the opposite direction of the light

% plot population vector (Fig. 2D)
subplot(2,2,4);
rurchin = 5; % size of the sea urchin
rmax = net.param.rPRCmax/2; % maximum r tick
h1 = polarplot([pop_vect(2)*pi/180,pop_vect(2)*pi/180],rurchin+[0, pop_vect(1)],'color',col(2,:),'linewidth',3);
x = rlim();
rmax = max([rmax,x(2),7]);
hold on;
h2 = polarplot(linspace(0,2*pi+0.1,500),rurchin+theta_p*ones(1,500),'color',col(1,:),'linewidth',2.5);
[ah] = ploturchin(gca, stim, rurchin, 4*rurchin, fntsz);

title('population vector');

hold off;

legend([h1,h2],'v_{pop}','\theta_p');
legend 'boxoff';

clear h1 h2;

sgtitle(sprintf('%i deg %s, Delta rho = %i deg, delta = %i deg',stim.width,stim_pattern,net.param.Delta_rho,net.param.delta), 'fontsize', fntsz+2);
