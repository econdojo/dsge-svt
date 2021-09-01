function nwhac = NeweyWest(H,nlag)
% Function NEWEYWEST
%
% Purpose:    Compute Newey-West heteroskedasticity and autocorrelation
%             consistent (HAC) covariance matrix estimator; see
%             Newey & West (1987,1994) for details
%
% Format:     nwhac = NeweyWest(H,nlag)
%
% Input:      H         demeaned data with ROW observations
%             nlag      number of lags
%
% Output:     nwhac     Newey-West HAC estimator
%
% Written by Fei Tan, Saint Louis University
% Updated: October 20, 2017

%% -------------------------------------------
%                HAC Estimator
%---------------------------------------------

% Initialization
M = size(H,1);                    % sample size
if nargin<2 || nlag<0
    nlag = floor(4*(M/100)^(2/9));% Newey & West (1994) automatic procedure
elseif nlag>floor(M/2)
    error('Number of lags exceeds half of sample size.')
end

% Compute HAC estimator
nwhac = H'*H/M;
for k = 1:nlag
    wk = 1-k/(nlag+1);            % Bartlett kernel
    Omega = H((k+1):end,:)'*H(1:(end-k),:)/M;
    nwhac = nwhac+wk*(Omega+Omega');
end

%-------------------- END --------------------