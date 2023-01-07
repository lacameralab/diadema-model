% This script generates Figure 3B in the paper.
%
% Tianshu Li
% Oct. 28, 2022

clear;

if ~exist('creatStimulus','file') addpath('./basic_functions/'); end

acc = 0.1; % degress, bin size for stimulus
stimwidth = 40; % deg, phi_stim
stimcenter = 180; % deg. Notice that the stimulus always centers at 0 deg in the simulation.

figure(1); clf;
subplot(6,1,1);
x = creatStimulus('bar', stimwidth, stimcenter, acc);
imagesc(x.intensity);  colormap 'gray';
colorbar 'off';
set(gca,'xtick',[]);
set(gca,'ytick',[]);
subplot(6,1,2);
x = creatStimulus('DoG', stimwidth, stimcenter, acc);
imagesc(x.intensity);  colormap 'gray';
colorbar 'off';
set(gca,'xtick',[]);
set(gca,'ytick',[]);
subplot(6,1,3);
x = creatStimulus('hermitian', stimwidth, stimcenter, acc);
imagesc(x.intensity);  colormap 'gray';
colorbar 'off';
set(gca,'xtick',[]);
set(gca,'ytick',[]);
subplot(6,1,4);
x = creatStimulus('square', stimwidth, stimcenter, acc);
imagesc(x.intensity);  colormap 'gray';
colorbar 'off';
set(gca,'xtick',[]);
set(gca,'ytick',[]);
subplot(6,1,5);
x = creatStimulus('2square', stimwidth, stimcenter, acc);
imagesc(x.intensity);  colormap 'gray';
colorbar 'off';
set(gca,'xtick',[]);
set(gca,'ytick',[]);
subplot(6,1,6);
x = creatStimulus('morlet', stimwidth, stimcenter, acc);
imagesc(x.intensity);  colormap 'gray';
colorbar 'off';
set(gca,'xtick',[]);
set(gca,'ytick',[]);