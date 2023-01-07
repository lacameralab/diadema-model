function [net,niter] = RNoutput(net,maxiter)
% [net,niter] = RNoutput(net,maxiter)
% This function calculates the output of all RNs given a network and
% PRC outputs.
% Inputs:
%   maxiter: maximum number of iterations. Default: 5000. 
% Outputs:
%   net: structure, contains network parameters
%   niter: number of iteration for recurrent connection to become stable
%
% Tianshu Li
% May 25, 2021

if nargin < 2
    maxiter = 5000; % maximum iterations for the recurrent net
end

if ~isfield(net.param,'beta_RN')
    beta_RN = 3; % Steepness of RN3s sigmoidal response function
    net.beta_RN = beta_RN;
else
    beta_RN = net.param.beta_RN;
end
if ~isfield(net.param,'xc_RN')
    xc_RN = -0.6; % Location parameter of RN's sigmoidal response function
    net.param.xc_RN = xc_RN;
else
    xc_RN = net.param.xc_RN;
end

% RN's output
eps = 1e-5;
for kAmb = 1:net.nAmbulacrum
    r_out = net.param.rRNmax*ones(net.RN(kAmb).nCell,1);
    for iter = 1:maxiter
        r_in_PRC = net.RN(kAmb).W_RP*net.PRC(kAmb).out; % feedforward input from PRC
        r_in_RN = net.RN(kAmb).W_RR*r_out; % recurrent input from RN
        r_in = r_in_PRC+r_in_RN; % total input
        r_out_next = net.RN(kAmb).rmax.*(1+tanh(beta_RN*(r_in./abs(net.RN(kAmb).W_RP*net.param.rPRCmax*ones(size(net.PRC(kAmb).out)))-xc_RN)))/2; % RN output, sec 5.3.3, eq. 7
        r_out_next(r_out_next<0) = 0; % cell activity cannot go below 0
        if norm(r_out_next-r_out) < eps % check convergence
            net.RN(kAmb).out = r_out_next;
            niter = iter;
            break;
        else
            if iter == maxiter
                warning('Ambulacrum %i does not converge!\n',kAmb);
            end
            r_out = r_out_next;
        end
    end
end

end