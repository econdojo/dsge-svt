% Function USER_SSP
%
% Purpose:    User input - steady state/implied parameters; see
%             Smets & Wouters (2007) for model details
%
% Format:     user_ssp
%
% Input:      P         parameter structure
%             para      model parameters
%
% Output:     ssp       steady state/implied parameters
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% User Input Start Here

% Steady state/implied parameters
ssp(P.delta) = 0.025;
ssp(P.lambda_w) = 1.50;
ssp(P.g_y) = 0.18;
ssp(P.eps_w) = 10;
ssp(P.eps_p) = 10;
ssp(P.beta) = 1/(para(P.beta_hat)/100+1);
ssp(P.gamma) = para(P.gamma_bar)/100+1;
ssp(P.pi_star) = para(P.pi_bar)/100+1;
ssp(P.r_bar) = 100*(ssp(P.gamma)^para(P.sigma_c)*ssp(P.pi_star)/ssp(P.beta)-1);
ssp(P.r_k_ss) = ssp(P.gamma)^para(P.sigma_c)/ssp(P.beta)-1+ssp(P.delta);
ssp(P.w_ss) = ((1-para(P.alpha))^(1-para(P.alpha))*para(P.alpha)^para(P.alpha)/para(P.Phi)/ssp(P.r_k_ss)^para(P.alpha))^(1/(1-para(P.alpha)));
ssp(P.i_k) = ssp(P.gamma)-1+ssp(P.delta);
ssp(P.l_k) = (1-para(P.alpha))/para(P.alpha)*ssp(P.r_k_ss)/ssp(P.w_ss);
ssp(P.k_y) = para(P.Phi)*ssp(P.l_k)^(para(P.alpha)-1);
ssp(P.i_y) = ssp(P.i_k)*ssp(P.k_y);
ssp(P.c_y) = 1-ssp(P.g_y)-ssp(P.i_y);
ssp(P.z_y) = ssp(P.r_k_ss)*ssp(P.k_y);
ssp(P.wl_c) = 1/ssp(P.lambda_w)*(1-para(P.alpha))/para(P.alpha)*ssp(P.r_k_ss)*ssp(P.k_y)/ssp(P.c_y);

%% User Input End Here