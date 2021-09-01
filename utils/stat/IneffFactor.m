function ineff = IneffFactor(chain,acf_nlag)
% Function INEFFFACTOR
%
% Purpose:    Compute inefficiency factor of a MCMC sample based on Parzen
%             window; see
%             Paolo Giordani, Michael Pitt, & Robert Kohn (2011): Bayesian
%             inference for time series state space models, The Oxford Handbook
%             of Bayesian Econometrics
%
% Format:     ineff = IneffFactor(chain,acf_nlag)
%
% Input:      chain     univariate time series (column vector)
%             acf_nlag  ACF truncation level; number of lags if >=1 integer
%
% Output:     ineff     inefficiency factor
%
% Written by Fei Tan, Saint Louis University
% Updated: July 18, 2017

%% -------------------------------------------
%             Inefficiency Factor
%---------------------------------------------

% Truncated autocorrelation function
M = length(chain);
acf = autocorr(chain,M-1);
if acf_nlag<1
    ind = find(abs(acf)<=acf_nlag,1);  % truncate small ACF values
    if isempty(ind)
        acf_nlag = floor(M/2);
    else
        acf_nlag = min([ind-1 floor(M/2)]);
    end
elseif acf_nlag>floor(M/2)
    error('Number of lags exceeds half of sample size.')
end
acf = acf(1:acf_nlag+1);

%Calculate Parzen weights
parzen = zeros(acf_nlag+1,1);
for k = 0:acf_nlag
    x = k/acf_nlag;
    if x<=0.5
        parzen(k+1) = 1-6*x^2+6*x^3;
    else
        parzen(k+1) = 2*(1-x)^3;
    end
end

% Calculate inefficiency factor
ineff = 2*sum(parzen.*acf)-1;

%-------------------- END --------------------
