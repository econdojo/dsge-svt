function [E_y,E_yy] = Moment(SSR,nlag)
% Function MOMENT
%
% Purpose:    Compute DSGE population moments; see
%             Del Negro & Schorfheide (2004), Priors from general equilibrium
%             models for VARs, International Economic Review
%
% Format:     [E_y,E_yy] = Moment(SSR,nlag)
%
% Input:      SSR       state space representation (structure)
%
%                       s(t) = C + G*s(t-1) + M*e(t), E[ee'] = Sigma_e (diagonal)
%                       y(t) = D + Z*s(t) + u(t), E[uu'] = Sigma_u (diagonal)
%
%             nlag      number of lags
%
% Output:     E_y       E[y(t)]
%             E_yy      E[y(t)*y(t-k)'] for k = 0,1,...,nlag
%
% Written by Fei Tan, Saint Louis University
% Updated: June 10, 2022

%% -------------------------------------------
%           DSGE Population Moments
%---------------------------------------------

% Initialization
[dimy,dims] = size(SSR.Z);        % numbers of observables & states
E_yy = zeros(dimy,dimy,nlag+1);   % profile of E[y(t)*y(t-k)']

% Invariant distribution
E_s = (eye(dims)-SSR.G)\SSR.C;
E_ss = dlyap_sym(SSR.G,SSR.M*diag(SSR.Sigma_e)*SSR.M',1);

% Compute moments
E_y = SSR.D+SSR.Z*E_s;
E_yy(:,:,1) = SSR.Z*E_ss*SSR.Z'+SSR.D*SSR.D'+diag(SSR.Sigma_u);
for k = 1:nlag
    E_yy(:,:,k+1) = SSR.Z*SSR.G^k*E_ss*SSR.Z'+SSR.D*SSR.D';
end

%-------------------- END --------------------