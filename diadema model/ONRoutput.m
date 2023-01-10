function [ONR,niter] = ONRoutput(net,maxiter)
% [ONR,niter] = ONRoutput(net,ONR,maxiter)
% This function calculates the output of all ONRs given a network and
% RN outputs.
% Inputs:
%   net: structure, contains network parameters
%   maxiter: maximum number of iterations. Default: 5000. 
% Outputs:
%   niter: number of iteration for recurrent connection to become stable
% 
% Tianshu Li
% May 25, 2021

if nargin < 3
    maxiter = 1e3; % maximum iterations for recurrent net
end

ONR = net.ONR;

if ~isfield(net.param,'beta_ONR')
    beta_ONR = 4.5; % Steepness of ONR's sigmoidal response function
    net.param.beta_ONR = beta_ONR;
else
    beta_ONR = net.param.beta_ONR;
end
if ~isfield(net.param,'xc_ONR')
    xc_ONR = -0.45; % Location parameter of ONR's sigmoidal response function
    net.param.xc_ONR = xc_ONR;
else
    xc_ONR = net.param.xc_ONR;
end

% put the output of RN in different ambulacrum together
r_out_RN = [];
for kAmb = 1:net.nAmbulacrum
    r_out_RN = [r_out_RN;net.RN(kAmb).out];
end

% ONR's output
eps = 1e-5;
r_out = net.param.rONRmax*ones(ONR.nCell,1); % output of ONR cells
for iter = 1:maxiter
    r_in_RN = ONR.W_OR*r_out_RN; % feedforward input from RN
    r_in_ONR = ONR.W_OO*r_out; % recurrent input from ONR
    r_in = r_in_RN+r_in_ONR; % total input to ONR
    r_out_next = ONR.rmax.*(1+tanh(beta_ONR*(r_in./abs(ONR.W_OR*net.param.rRNmax*ones(size(r_in_RN)))-xc_ONR)))/2; % output, sec. 5.3.4, eq. 10
    r_out_next(r_out_next<0) = 0; % cell activity cannot go below 0
    if norm(r_out_next-r_out) < eps % check convergence
        ONR.out = r_out_next;
        niter = iter;
        break;
    else
        if iter == maxiter
            error('ONR does not converge!\n');
        end
        r_out = r_out_next;
    end
end

end