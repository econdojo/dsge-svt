function interval = ProbBand(chain,prob,hpdi)
% Function POSTSTAT
%
% Purpose:    Compute posterior probability bands; see
%             Chen & Shao (1999), Monte carlo estimation of Bayesian credible
%             and HPD intervals, J of Computational & Graphical Statistics
%
% Format:     interval = ProbBand(chain,prob,hpdi)
%
% Input:      chain     MCMC posterior draws
%             prob      credible set mass
%             hpdi      highest prob density or equal-tail interval
%
% Output:     interval  probability bands
%
% Written by Fei Tan, Saint Louis University
% Updated: June 15, 2017

%% -------------------------------------------
%         Posterior Probability Bands
%---------------------------------------------

% Initialization
[M,npara] = size(chain);     % sample size & number of parameters
order = sort(chain);         % sorted chain

% Compute probability bands
if hpdi                      % highest prob density interval (assuming unimodal)
    nin = round(M*prob);
    nout = M-nin+1;
    dist = order(nin:M,:)-order(1:nout,:);
    [dist,ind] = min(dist);
    ind = sub2ind(size(order),ind,1:npara);
    interval = [order(ind);order(ind)+dist]';
else                         % equal-tail interval
    interval = order([round(M*(1-prob)/2) round(M-M*(1-prob)/2)],:)';
end

%-------------------- END --------------------
