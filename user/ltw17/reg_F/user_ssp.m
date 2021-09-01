% Function USER_SSP
%
% Purpose:    User input - steady state/implied parameters
%
% Format:     user_ssp
%
% Input:      P         parameter structure
%             para      model parameters
%
% Output:     ssp       steady state/implied parameters
%
% Written by Fei Tan, Saint Louis University
% Updated: August 24, 2017

%% User Input Start Here

% Steady state parameters
% Calibration
ssp(P.mat) = 20;
ssp(P.beta) = 0.99;
ssp(P.delta) = 0.025;
ssp(P.alpha) = 0.33;
ssp(P.eta_p) = 0.14;
ssp(P.eta_w) = 0.14;
ssp(P.sgy) = 0.11;      % fiscal variables: U.S. average over 1954:Q3-2014:Q2
ssp(P.sby) = 1.47;
ssp(P.tau_c) = 0.023;
ssp(P.tau_k) = 0.218;
ssp(P.tau_l) = 0.186;
% Fixed
ssp(P.mu_h) = 0;
ssp(P.gam_g) = 0;
ssp(P.gam_k) = 0;
ssp(P.gam_l) = 0;
ssp(P.gam_z) = 0;
ssp(P.rho_k) = 0;
ssp(P.rho_l) = 0;
ssp(P.rho_c) = 0;
% Steady state
ssp(P.rho) = (1-1/ssp(P.mat))*(1/ssp(P.beta));
ssp(P.gam) = para(P.gam100)/100;
ssp(P.expg) = exp(ssp(P.gam));
ssp(P.pi_ss) = para(P.pbar)/100+1;
ssp(P.R_ss) = ssp(P.expg)*ssp(P.pi_ss)/ssp(P.beta);
ssp(P.R_bar) = (ssp(P.R_ss)/ssp(P.pi_ss)-1)*100;
ssp(P.p_b) = 1/(ssp(P.R_ss)-ssp(P.rho));
ssp(P.R_k) = (ssp(P.expg)/ssp(P.beta)-1+ssp(P.delta))/(1-ssp(P.tau_k));
ssp(P.r_k) = ssp(P.R_k);
ssp(P.psi1) = ssp(P.R_k)*(1-ssp(P.tau_k));
ssp(P.mc) = 1/(1+ssp(P.eta_p));
ssp(P.w) = (ssp(P.mc)*((1-ssp(P.alpha))^(1-ssp(P.alpha)))*(ssp(P.alpha)^ssp(P.alpha))*(ssp(P.R_k)^(-ssp(P.alpha))))^(1/(1-ssp(P.alpha)));
ssp(P.kl) = (ssp(P.w)/ssp(P.R_k))*ssp(P.alpha)/(1-ssp(P.alpha));
ssp(P.omegl) = (ssp(P.kl)^ssp(P.alpha))-ssp(P.R_k)*(ssp(P.kl))-ssp(P.w);
ssp(P.yl) = (ssp(P.kl)^ssp(P.alpha))-ssp(P.omegl);
ssp(P.il) = (1-(1-ssp(P.delta))*exp(-ssp(P.gam)))*ssp(P.expg)*ssp(P.kl);
ssp(P.cl) = ssp(P.yl)*(1-ssp(P.sgy))-ssp(P.il);
ssp(P.zl) = ((1-ssp(P.R_ss)/ssp(P.pi_ss)*exp(-ssp(P.gam)))*ssp(P.sby)-ssp(P.sgy))*ssp(P.yl)+ssp(P.tau_c)*ssp(P.cl)+ssp(P.tau_l)*ssp(P.w)+ssp(P.tau_k)*ssp(P.r_k)*ssp(P.kl);
ssp(P.c_nl) = ((1-ssp(P.tau_l))*ssp(P.w)+ssp(P.zl))/(1+ssp(P.tau_c));
ssp(P.c_sl) = (ssp(P.cl)-ssp(P.mu_h)*ssp(P.c_nl))/(1-ssp(P.mu_h));
ssp(P.c_starl) = ssp(P.c_sl) + para(P.alpha_g)*ssp(P.sgy)*ssp(P.yl);
ssp(P.l) = ((ssp(P.w)*(1-ssp(P.tau_l))/((1+ssp(P.tau_c))*(1+ssp(P.eta_w))))*(1/((1-para(P.theta)*exp(-ssp(P.gam)))*ssp(P.c_starl))))^(1/(para(P.xi)+1));
ssp(P.c_s) = ssp(P.c_sl)*ssp(P.l);
ssp(P.c_n) = ssp(P.c_nl)*ssp(P.l);
ssp(P.y) = ssp(P.yl)*ssp(P.l);
ssp(P.k) = ssp(P.kl)*ssp(P.l);
ssp(P.omeg) = ssp(P.omegl)*ssp(P.l);
ssp(P.c) = ssp(P.cl)*ssp(P.l);
ssp(P.inv) = ssp(P.il)*ssp(P.l);
ssp(P.z) = ssp(P.zl)*ssp(P.l);
ssp(P.b) = ssp(P.sby)*ssp(P.y);
ssp(P.g) = ssp(P.sgy)*ssp(P.y);
ssp(P.ky) = ssp(P.k)/ssp(P.y);
ssp(P.cy) = ssp(P.c)/ssp(P.y);
ssp(P.ly) = ssp(P.l)/ssp(P.y);
ssp(P.S) = ssp(P.tau_k)*ssp(P.r_k)*ssp(P.k)+ssp(P.tau_l)*ssp(P.l)*ssp(P.w)+ssp(P.tau_c)*ssp(P.c)-ssp(P.g)-ssp(P.z);
ssp(P.lamp) = ((1+ssp(P.beta)*para(P.chi_p))*para(P.omega_p))/((1-ssp(P.beta)*para(P.omega_p))*(1-para(P.omega_p)));
ssp(P.lamw) = para(P.omega_w)*(1+ssp(P.beta))*(1+para(P.xi)*(1+1/ssp(P.eta_w)))/((1-para(P.omega_w)*ssp(P.beta))*(1-para(P.omega_w)));

%% User Input End Here