function [Y,X] = DatMat(data,nlag)
% Function DATMAT
%
% Purpose:    Construct data matrices for vector autoregression (VAR)
%
% Format:     [Y,X] = DatMat(data,nlag)
%
% Input:      data      data set with ROW observations
%             nlag      number of VAR lags
%
% Output:     Y         data matrix of dependent variables
%             X         data matrix of independent variables
%
% Written by Fei Tan, Saint Louis University
% Updated: July 18, 2017

%% -------------------------------------------
%           Construct Data Matrices
%---------------------------------------------

Y = data((nlag+1):end,:);    % exclude initial p lags
T = size(Y,1);               % time span
X = ones(T,1);               % regressors
for k = 1:nlag
    X = cat(2,X,data((nlag-k+1):(nlag-k+T),:));
end

%-------------------- END --------------------