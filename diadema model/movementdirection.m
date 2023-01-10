function [movedirect] = movementdirection(v_pop,theta_p,predirection)
% [movedirect] = movementdirection(v_pop,theta_p,predirection,sigma_vpop)
% This function estimates the movement direction of the animal given the
% population vector v_pop.
% The new movement direction is sampled from a mixture of 2 Normal
% distribution one centered at the direction of v_pop, the other centered
% at the previous movement direction (phi_pre).
% See section 5.3.6 for detail.
% Inputs:
%   v_pop: population vector, v_pop = [angle, length].
%   theta_p: threshold for visibility
%   predirection: phi_pre, the movement direction of the previous step, unit:
%       degree. When there is no previous step, predirection is NaN, and
%       the corresponding ditribution is uniform distribution.
% Output:
%   movedirect: movement direction, unit: degree
% 
% Tianshu Li
% Nov. 11, 2021

f = @(x) 1./(1+exp(-10*x)); % sigmoid function, mixing proportion q = f(|v_pop|-theta_p) (Eq. 12)

q = f(v_pop(:,2)-theta_p); % mixing proportion
mu1 = v_pop(:,1);
epsilon = 1e-5; % to keep sigma1 positive
sigma1 = max(1./((v_pop(:,2)-theta_p)+epsilon),epsilon);
mu2 = predirection;
sigma2 = 10*ones(size(sigma1)); % fixed sigma_2

x = rand(size(q)); % draw x from categorical distribution parametrized by vector q
idx1 = find(x < q); % the ones follow the first Normal distribution
idx2 = find(x >= q); % the ones follow the second Normal distribution
movedirect = zeros(length(q),1);
movedirect(idx1) = mu1(idx1)+sigma1(idx1).*randn(length(idx1),1); % draw samples from the first Normal distribution
if ~isnan(mu2)
    movedirect(idx2) = mu2+sigma2(idx2).*randn(length(idx2),1); % draw samples from the second Normal distribution
else
    movedirect(idx2) = 360*rand(length(idx2),1); % first step (no previous step), draw samples from the second distribution (uniform distribution)
end

movedirect = mod(movedirect,360);

end