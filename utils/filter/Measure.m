function y_prob = Measure(y,SSR,swarm,pred)
% Function MEASURE
%
% Purpose:    Compute likelihoods of observables from measurement equation
%
% Format:     y_prob = Measure(y,SSR,swarm,pred)
%
% Input:      y         current observation (row vector)
%             SSR       state space representation (structure)
%
%                       y(t) = Measure(s(t)) + u(t), E[uu'] = Sigma_u (diagonal)
%
%             swarm     collection of particles
%             pred      current swarm = 0; previous swarm = 1
%
% Output:     y_prob    log likelihoods of current observation
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%           Measurement Likelihoods            
%---------------------------------------------

% Initialization
N = size(swarm,2);           % numbers of particles
dimy = length(y);            % number of observables
c = 1;                       % scaling parameter
pdof = 15;                   % proposal t degrees of freedom

% Compute likelihoods of current observation
if pred
    swarm = repmat(SSR.C,1,N)+SSR.G*swarm;
    Sigma_u = c*(SSR.Z*SSR.M*diag(SSR.Sigma_e)*SSR.M'*SSR.Z'+diag(SSR.Sigma_u));
    y_prob = mvt_pdf(repmat(y,N,1)-repmat(SSR.D',N,1)-swarm'*SSR.Z',zeros(1,dimy),Sigma_u,pdof);
else
    y_prob = mvt_pdf(repmat(y,N,1)-repmat(SSR.D',N,1)-swarm'*SSR.Z',zeros(1,dimy),diag(SSR.Sigma_u),inf);
end

%-------------------- END --------------------
