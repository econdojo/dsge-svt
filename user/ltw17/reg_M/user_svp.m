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
% svp_ss(V.e_a) = P.sigma_a;
% svp_ss(V.e_b) = P.sigma_b;
% svp_ss(V.e_i) = P.sigma_i;
% svp_ss(V.e_p) = P.sigma_p;
% svp_ss(V.e_w) = P.sigma_w;
% svp_ss(V.e_m) = P.sigma_m;
% svp_ss(V.e_g) = P.sigma_g;
% svp_ss(V.e_z) = P.sigma_z;
% stochastic volatility
svp_ss(V.e_a) = P.mu_a;
svp_ss(V.e_b) = P.mu_b;
svp_ss(V.e_i) = P.mu_i;
svp_ss(V.e_p) = P.mu_p;
svp_ss(V.e_w) = P.mu_w;
svp_ss(V.e_m) = P.mu_m;
svp_ss(V.e_g) = P.mu_g;
svp_ss(V.e_z) = P.mu_z;

% ar indices
svp_ar(V.e_a) = P.phiv_a;
svp_ar(V.e_b) = P.phiv_b;
svp_ar(V.e_i) = P.phiv_i;
svp_ar(V.e_p) = P.phiv_p;
svp_ar(V.e_w) = P.phiv_w;
svp_ar(V.e_m) = P.phiv_m;
svp_ar(V.e_g) = P.phiv_g;
svp_ar(V.e_z) = P.phiv_z;

% vv indices
svp_vv(V.e_a) = P.sig2_a;
svp_vv(V.e_b) = P.sig2_b;
svp_vv(V.e_i) = P.sig2_i;
svp_vv(V.e_p) = P.sig2_p;
svp_vv(V.e_w) = P.sig2_w;
svp_vv(V.e_m) = P.sig2_m;
svp_vv(V.e_g) = P.sig2_g;
svp_vv(V.e_z) = P.sig2_z;

%% User Input End Here