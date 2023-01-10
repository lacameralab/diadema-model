% This script simulates the trajectories (Fig. 5).
% See section 5.3.6 for detail.
%
% P.S.: the unit of angles in this script is "degree" ([0, 360]).
%
% Tianshu Li
% June 17, 2021

clear;

addpath('./utils/');

%% set stimulus
theta_p = 5; % threshold for |v_pop| (length of population vector)
Ntrials = 100; % number of trials
dx = 0.1; % step size
maxnstep = 100; % maximum number of steps in one trial, if reaches the number but the animal still have not reach the wall, force the trial to stop
stop_radius = 3/4; % when the center of the animal is this value away from the center of the arena, the animal stops (1-stop_radius is the radius of the animal).

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

%% trajectory (Ntrials animals in arena, random initial orientation)
stimcenters = 72*rand([1,Ntrials]); % pentaradial symmetry, stimulus center between 0 and 72 is the same as between 0 and 360
animal_trajectory = {}; % each entry is a matrix with the positions of the animal at each step in one trial

for ktrial = 1:Ntrials
    % generate stimulus with respect to animal's initial orientation
    stim.center = stimcenters(ktrial); % deg, center of stimulus
    if strcmp(stim_pattern,'control')
        stim.intensity = 0.77*ones(size(stim.intensity)); % control
        stim0 = stim;
        stim0.intensity = 0.77*ones(size(stim.intensity)); % control
    else
        stim = createStimulus(stim_pattern,stim.width,stim.center,acc);
        stim0 = createStimulus(stim_pattern,stim.width,0,acc);
    end

    animal_orientation = -stim.center; % the direction where the first ambulacrum faces in arena s coordinate system (Fig. S1)
    animalposition = [animal_orientation,0]; %(theta, r), the polar coordinates of the current position of animal (scale the radius of the arena to be 1).
    animalposition_next = animalposition; %(theta, r), the polar coordinates of the next position of animal
    animal_trajectory{ktrial} = [];
    dmovement = [0,dx]; % (theta, r), the movement arrow each step (vector delta r)

    nstep = 0; % number of steps
    while animalposition_next(2) < stop_radius && nstep < maxnstep % the animal has not reach the wall or the maximum number of steps
        animalposition = animalposition_next;
        animal_trajectory{ktrial} = [animal_trajectory{ktrial};animalposition]; % add the current position to trajectory

        newstim = stimatposition(stim0,stim,animalposition);
        net = PRCoutput(net,newstim);
        net = RNoutput(net);
        net.ONR = ONRoutput(net);

        z = sum((net.ONR.out).*(cosd(net.ONR.phi_pref) + sind(net.ONR.phi_pref)*1i))/sqrt(net.ONR.nCell); % estimated direction of light by the animal
        z_amplitude = abs(z);
        z_direction = mod(animal_orientation+rad2deg(angle(z)),360); % movement direction

        % movement direction and the position of the animal after this step
        if nstep == 0
            if z_amplitude<=theta_p % if this is the first step, move randomly.
                % uniform distribution if the amplitude of population vector is smaller than threshold theta_p
                dmovement(1) = 360*rand();
            else
                dmovement(1) = movementdirection([z_direction,z_amplitude],theta_p,NaN);
            end
        else
            dmovement(1) = movementdirection([z_direction,z_amplitude],theta_p,dmovement(1));
        end

        animalposition_next_complex = dmovement(2)*(cosd(dmovement(1))+1i*sind(dmovement(1))) + animalposition(2)*(cosd(animalposition(1))+1i*sind(animalposition(1)));
        animalposition_next(1) = mod(rad2deg(angle(animalposition_next_complex)),360);
        animalposition_next(2) = abs(animalposition_next_complex);

        nstep = nstep +1;

        if animalposition_next(2) > stop_radius % last step, this step ends when the animal hits the wall
            animalposition_complex = animalposition(2)*(cosd(animalposition(1))+1i*sind(animalposition(1)));
            dmovement_complex = dmovement(2)*(cosd(dmovement(1))+1i*sind(dmovement(1)));
            dmovement_unit_complex = dmovement_complex/abs(dmovement_complex); % the unit vector of the movement
            x = real(animalposition_complex);
            y = imag(animalposition_complex);
            e1 = real(dmovement_unit_complex);
            e2 = imag(dmovement_unit_complex);
            k = (-2*(x*e1+y*e2)+2*sqrt((x*e1+y*e2)^2-(e1^2+e2^2)*(x^2+y^2-stop_radius^2)))/(2*(e1^2+e2^2)); % length of the last step
            animalposition_next_complex = animalposition_complex+k*dmovement_unit_complex;
            animalposition_next(1) = mod(rad2deg(angle(animalposition_next_complex)),360);
            animalposition_next(2) = abs(animalposition_next_complex)+1e-5; % add a small value 1e-5 to avoid numerical issues which leads to animalposition_next(2)<stop_radius

            animalposition = animalposition_next;
            animal_trajectory{ktrial} = [animal_trajectory{ktrial}; animalposition]; % add the current position to trajectory
        end
    end

end

%% visualize the trajectory
figure(1); clf;
set(gcf,'position',[89    66   845   649]);
stim = createStimulus(stim_pattern,stim.width,0,acc);
if strcmp(stim_pattern,'control')
    stim.intensity  = 0.77*ones(size(stim.intensity)); % control
end

arenawall = 1; % radius of the arena wall
for i = 1:10
    [strike,dip]= meshgrid(stim.phi*pi/180,arenawall+(i-1)*0.01);
    polarscatter(strike,dip,20*ones(size(strike)),repmat(stim.intensity',1,3),'filled','MarkerFaceAlpha',0.8);
    hold on;
end
polarplot(stim.phi*pi/180,(arenawall-0.01)*ones(size(stim.phi)),'k','linewidth',2);
polarplot(stim.phi*pi/180,(arenawall+0.06)*ones(size(stim.phi)),'k','linewidth',2);
if arenawall > stop_radius
    polarplot(stim.phi*pi/180,stop_radius*ones(size(stim.phi)),'k--','linewidth',2);
end
colormap bone;

pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'counterclockwise';
pax.FontSize = 22;
rticklabels('');
rlim([0,arenawall+0.06]);

hold on;
% plot trajectory
ntrial_highlight = ceil(Ntrials/10); % choose 10% of the trials to highlight
trial_highlight_idx = randperm(Ntrials,ntrial_highlight); % put the trial indices which you want to hightlight here
for ktrial = 1:Ntrials
    h2 = polarplot(deg2rad(animal_trajectory{ktrial}(:,1)),animal_trajectory{ktrial}(:,2),'linewidth',2);

    if ~any(trial_highlight_idx==ktrial)
        h2.Color = [h2.Color,0.1];
    end

    h3 = polarscatter(deg2rad(animal_trajectory{ktrial}(end,1)),animal_trajectory{ktrial}(end,2),100,h2.Color,'filled');
    h4 = polarscatter(deg2rad(animal_trajectory{ktrial}(end,1)),arenawall-0.01,100,0.6*[1,1,1],'filled');
    h4.MarkerEdgeColor = [0,0,0];
end

title(sprintf('%i trials, %i %s, Delta rho = %i deg, delta = %i deg',Ntrials,stim.width,stim_pattern,net.param.Delta_rho,net.param.delta),'fontsize',20);
