function log_vol = PostVol(theta,data,mixid)
% Function POSTVOL
%
% Purpose:    Draw posterior volatilities
%
% Format:     log_vol = PostVol(theta,data,mixid)
%
% Input:      theta     sv parameters (current)
%             data      transformed data (current)
%             mixid     mixture normal indices (current)
%
% Output:     log_vol   log volatilities (next)
%
% Written by Fei Tan, Saint Louis University
% Updated: January 15, 2020

%% -------------------------------------------
%             Sample Volatilities
%---------------------------------------------

% Initialization
[p10,m10,v10] = mix10;       % 10-component mixture normal
mu = theta(1);               % sv parameters
phi = theta(2);
sig2 = theta(3);

% Sample log volatilities
for k = 1:10
    SSR.C{k} = (1-phi)*mu;
    SSR.G{k} = phi;
    SSR.Sigma_e{k} = sig2;
    SSR.D{k} = m10(k);
    SSR.Z{k} = 1;
    SSR.Sigma_u{k} = v10(k);
end
[~,~,~,log_vol,~] = GibbsSmoother_mex(data,SSR,repmat(p10,10,1),mixid,true,false);

%-------------------- END --------------------