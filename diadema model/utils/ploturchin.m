function [h] = ploturchin(gca, stim, urchinboundary, arenaboundary, fntsz, psi)
% [h,ah] = ploturchin(gca, stim, urchinboundary, arenaboundary, stimwidth, fntsz, psi)
% Inputs:
%   stim: stimulus, structure containing phi (angle) and intensity (light
%   intensity corresponding to each angle in phi).
%   urchinboundary: r-value of the boundary of the urchin
%   arenaboundary: r-value of the boundary of arena
%   fntsz: fontsize
%   stimwidth: width of stimulus
%   psi: orientation of the animal (direction of the first ambulacrum in
%       the arena). unit: degree.
% Output:
%   h: line and scatter handle

if nargin < 5
    fntsz = 18;
end

if nargin < 6
    psi = 0;
end

% black circle represent the boundary of the animal
h(1) = polarplot(stim.phi*pi/180,urchinboundary*ones(size(stim.intensity)),'k','linewidth',2);
hold on;

% ambulacra
for kAmb = 1:5
    h(kAmb+1) = polarplot((2*pi/5*(kAmb-1)+psi*pi/180)*[1,1],[0,urchinboundary],'k','linewidth',3);
end

% stimulus
for i = 1:10
    [strike,dip]= meshgrid(stim.phi*pi/180,(9+(i-1)*0.1)*arenaboundary/10);
    polarscatter(strike,dip,50*ones(size(strike)),repmat(stim.intensity',1,3),'filled','MarkerFaceAlpha',0.8);
end
polarplot(stim.phi*pi/180,(9-0.2)*arenaboundary*ones(size(stim.phi))/10,'k','linewidth',2);
polarplot(stim.phi*pi/180,(10)*arenaboundary*ones(size(stim.phi))/10,'k','linewidth',2);
colormap bone;

rticks([urchinboundary (urchinboundary+arenaboundary)/2]);
rticklabels({'0',sprintf('%g',(arenaboundary-urchinboundary)/2)});

pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'counterclockwise';
pax.FontSize = fntsz;

end