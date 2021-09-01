% Function USER_SVP
%
% Purpose:    User input - stochastic volatility parameters
%
% Format:     user_svp
%
% Input:      P         parameter structure
%             V         variable structure
%
% Output:     svp_ss    sv steady state indices
%             svp_ar    AR parameter indices
%             svp_vv    vol-vol parameter indices
%
% Written by Fei Tan, Saint Louis University
% Updated: January 15, 2020

%% User Input Start Here

% ss indices
% constant volatility
% svp_ss(V.eps_R) = P.sigma_R;
% svp_ss(V.eps_G) = P.sigma_G;
% svp_ss(V.eps_Z) = P.sigma_Z;
% stochastic volatility
svp_ss(V.eps_R) = P.mu_R;
svp_ss(V.eps_G) = P.mu_G;
svp_ss(V.eps_Z) = P.mu_Z;

% ar indices
svp_ar(V.eps_R) = P.phi_R;
svp_ar(V.eps_G) = P.phi_G;
svp_ar(V.eps_Z) = P.phi_Z;

% vv indices
svp_vv(V.eps_R) = P.sig2_R;
svp_vv(V.eps_G) = P.sig2_G;
svp_vv(V.eps_Z) = P.sig2_Z;

%% User Input End Here